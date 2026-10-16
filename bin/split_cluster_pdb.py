#!/usr/bin/env python3
"""
bin/split_cluster_pdb.py
Splits GROMACS multi-model cluster PDB into centroid structures and in-situ ligand poses.
Extracts both PDB and MAE for centroids, and SDF and MAE for ligand poses.
Outputs structured receipt JSON for CLI DB synchronization.
"""

import argparse
import sys
from pathlib import Path
from schrodinger import structure


STANDARD_PROTEIN = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "CYX",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASH",
    "GLH",
    "ACE",
    "NME",
    "NMA",
}


def normalize_gmx_cluster_pdb(raw_pdb_path: Path, lig_name: str) -> list[Path]:
    """
    Parses multi-model GROMACS cluster PDB, resolves chain boundaries,
    injects canonical Chain IDs ('A', 'B'... and 'Z' for ligands), inserts strict TER cards,
    and emits clean, un-interleaved model PDBs.
    """
    lines = raw_pdb_path.read_text(errors="ignore").splitlines()
    models = []
    current_model_lines = []
    in_model = False

    for line in lines:
        if line.startswith("MODEL"):
            in_model = True
            current_model_lines = []
            continue
        elif line.startswith("ENDMDL"):
            in_model = False
            if current_model_lines:
                models.append(current_model_lines)
                current_model_lines = []
            continue
        elif in_model or line.startswith(("ATOM", "HETATM")):
            current_model_lines.append(line)

    if not models and current_model_lines:
        models.append(current_model_lines)

    clean_model_files = []
    for m_idx, m_lines in enumerate(models, 1):
        clean_lines = []
        chain_chars = "ABCDEFGHIJKLMNOPQRSTUVWXY"
        chain_idx = 0
        curr_chain = chain_chars[chain_idx]
        prev_res_seq = None
        prev_res_name = None
        global_atom_serial = 1

        for line in m_lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_name = line[12:16].strip()
            res_name = line[17:20].strip().upper()
            try:
                res_seq = int(line[22:26].strip())
            except ValueError as e:
                print(f"FATAL: Malformed residue sequence number in line '{line}': {e}", file=sys.stderr)
                sys.exit(1)

            # Chain transition detection rules
            is_lig = (res_name == lig_name) or (res_name not in STANDARD_PROTEIN)
            prev_is_lig = (prev_res_name == lig_name) or (
                prev_res_name not in STANDARD_PROTEIN if prev_res_name else False
            )

            is_chain_break = False
            if prev_res_seq is not None:
                # 1. Monotonic reset to 1 in GROMACS
                if res_seq < prev_res_seq:
                    is_chain_break = True
                # 2. Capped terminal transition NME -> ACE
                elif prev_res_name in ["NME", "NMA"] and res_name in ["ACE"]:
                    is_chain_break = True
                # 3. Transition to/from small molecule
                elif is_lig != prev_is_lig:
                    is_chain_break = True

            if is_chain_break:
                clean_lines.append(
                    f"TER   {global_atom_serial:5d}      {prev_res_name:>3s} {curr_chain}{prev_res_seq:4d}"
                )
                global_atom_serial += 1
                if is_lig:
                    curr_chain = "Z"
                else:
                    chain_idx += 1
                    curr_chain = chain_chars[chain_idx % len(chain_chars)]

            atom_field = (
                f"{atom_name:<4s}" if len(atom_name) == 4 else f" {atom_name:<3s}"
            )
            coords_and_rest = line[30:54]
            element = atom_name[:1] if not atom_name[:1].isdigit() else atom_name[1:2]

            rec_type = "HETATM" if is_lig else "ATOM  "
            new_line = (
                f"{rec_type}{global_atom_serial:5d} {atom_field} {res_name:>3s} "
                f"{curr_chain}{res_seq:4d}    {coords_and_rest:24s}  1.00  0.00          {element:>2s}"
            )
            clean_lines.append(new_line)
            global_atom_serial += 1
            prev_res_seq = res_seq
            prev_res_name = res_name

        if clean_lines:
            clean_lines.append(
                f"TER   {global_atom_serial:5d}      {prev_res_name:>3s} {curr_chain}{prev_res_seq:4d}"
            )
            clean_lines.append("END")
            out_file = raw_pdb_path.parent / f"norm_model_{m_idx}.pdb"
            out_file.write_text("\n".join(clean_lines) + "\n")
            clean_model_files.append(out_file)

    return clean_model_files


def main():
    parser = argparse.ArgumentParser(
        description="Split cluster PDB and extract ligand poses with rigorous topology & chain preservation"
    )
    parser.add_argument("--input", required=True, help="Input clusters_multi.pdb")
    parser.add_argument("--sys_name", required=True, help="Parent system name")
    parser.add_argument("--ligand", default="", help="Ligand residue name (e.g. MK1)")
    parser.add_argument(
        "--max_clusters", type=int, default=0, help="Max clusters to extract"
    )
    parser.add_argument(
        "--out_receipt", default="cluster_receipt.json", help="Output receipt JSON"
    )
    args = parser.parse_args()

    in_path = Path(args.input)
    if not in_path.exists() or in_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Cluster PDB {in_path} not found or empty.")

    import json

    lig_name = args.ligand.strip().upper()
    max_c = args.max_clusters
    receipt_records = []

    # 1. Pre-process GROMACS cluster PDB to inject canonical Chain IDs ('A', 'B'..., 'Z') and TER separators
    clean_model_pdbs = normalize_gmx_cluster_pdb(in_path, lig_name)

    from schrodinger.structutils import assignbondorders

    water_residues = {"HOH", "WAT", "SOL", "TIP", "TIP3", "TIP4", "W"}

    for i, model_pdb in enumerate(clean_model_pdbs, 1):
        if max_c > 0 and i > max_c:
            break

        st = structure.StructureReader.read(str(model_pdb))
        centroid_poses = {}

        # 2. Extract ALL non-protein, non-water small molecules (ligand and cofactors)
        small_mol_resnames = set()
        for res in st.residue:
            r_clean = res.pdbres.strip().upper()
            if r_clean not in STANDARD_PROTEIN and r_clean not in water_residues:
                small_mol_resnames.add(r_clean)

        c_title = f"{args.sys_name}_centroid_{i}"
        view_pdb_name = f"{c_title}_view.pdb"

        # 3. Clean AMBER residue translations across all residues
        res_translation = {
            "HIE": "HIS",
            "HID": "HIS",
            "HIP": "HIS",
            "CYX": "CYS",
            "CYM": "CYS",
            "ASH": "ASP",
            "GLH": "GLU",
            "LYN": "LYS",
        }
        for res in st.residue:
            r_clean = res.pdbres.strip().upper()
            if r_clean in res_translation:
                res.pdbres = res_translation[r_clean]

        # Save complete assembled view PDB of the centroid system for PyMOL
        st.title(c_title)
        st.property["s_m_title"] = c_title
        st.write(view_pdb_name)

        # 4. Extract ALL non-protein, non-water small molecules into separate poses
        for r_code in sorted(list(small_mol_resnames)):
            mol_atoms = [
                a.index
                for a in st.atom
                if a.pdbres.strip().upper() == r_code
            ]
            if mol_atoms:
                mol_st = st.extract(mol_atoms)
                
                if len(mol_st.molecule) > 1 and len(mol_atoms) > len(mol_st.molecule):
                    raise ValueError(
                        f"FATAL: Cluster centroid residue '{r_code}' resolves to {len(mol_st.molecule)} disconnected molecules! "
                        f"This corrupts topology. Ensure input complex had uniquely named molecules."
                    )
                
                mol_st.title = r_code
                mol_st.property["s_m_title"] = r_code

                # Perceive chemistry, double bonds and formal charges from 3D geometry
                assignbondorders.assign_st(mol_st)

                sdf_out = f"{r_code}_pose.sdf"
                mae_out = f"{r_code}_pose.mae"
                mol_st.write(sdf_out)
                mol_st.write(mae_out)
                st.deleteAtoms(mol_atoms)
                centroid_poses[r_code] = {
                    "sdf": sdf_out,
                    "mae": mae_out,
                }

        pdb_name = f"{c_title}_receptor.pdb"
        mae_name = f"{c_title}_receptor.mae"
        st.write(pdb_name)
        st.write(mae_name)

        receipt_records.append({
            "sys_name": args.sys_name,
            "centroid_title": c_title,
            "view_pdb": view_pdb_name,
            "pdb": pdb_name,
            "mae": mae_name,
            "poses": centroid_poses,
        })

        # Cleanup temporary normalized model PDB
        model_pdb.unlink(missing_ok=True)

    # Always write deterministic JSON receipt
    out_receipt_p = Path(args.out_receipt)
    json_path = out_receipt_p.with_suffix(".json") if out_receipt_p.suffix.lower() != ".json" else out_receipt_p
    json_path.write_text(json.dumps(receipt_records, indent=2), encoding="utf-8")

    print(f"Extracted {len(receipt_records)} cluster centroid(s) -> {json_path}")


if __name__ == "__main__":
    main()
