#!/usr/bin/env python3
"""
bin/extract_centroids.py
Atomic Stage 3: Direct seek frame extraction of physical medoid structures.
Preserves true chemical graph topology (formal charges, bond orders, stereochemistry)
via coordinate transfer onto pristine SubstanceState templates.
Extracts pure receptor, complex view, and all individual cofactor/ion poses.
"""

import argparse
import io
import json
import sys
from pathlib import Path
import MDAnalysis as mda
from rdkit import Chem
import numpy as np


# Inject framework root into sys.path to access canonical chem engine
FRAMEWORK_DIR = Path(__file__).resolve().parent.parent
if str(FRAMEWORK_DIR) not in sys.path:
    sys.path.insert(0, str(FRAMEWORK_DIR))

from cli.core.chem import transfer_coordinates_to_template


class NonClosingStringIO(io.StringIO):
    """
    In-memory string stream that preserves its buffer when closed by context managers.
    Allows MDAnalysis PDBWriter to cleanly finalize CONECT/END records without destroying the buffer.
    """
    def close(self):
        self.flush()


STANDARD_PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY", "HIS",
    "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL", "ASH", "GLH", "ACE", "NME", "NMA",
}

WATER_RESIDUES = {"HOH", "WAT", "SOL", "TIP", "TIP3", "TIP4", "W"}


def main():
    parser = argparse.ArgumentParser(description="Extract validated 3D medoid structures from trajectory")
    parser.add_argument("--struct", required=True, help="Clean reference structure PDB")
    parser.add_argument("--xtc", nargs="+", required=True, help="Clean trajectory XTC file(s)")
    parser.add_argument("--medoids_json", required=True, help="Input medoids_summary.json deliverable")
    parser.add_argument("--heterogens_spec", default="[]", help="JSON string with heterogens specification")
    parser.add_argument("--templates_dir", default="", help="Directory with staged template SDF files")
    parser.add_argument("--ligand_resname", default="", help="Active ligand residue name")
    parser.add_argument("--outdir", default="centroids_out", help="Output directory for structures")
    args = parser.parse_args()

    struct_p = Path(args.struct).resolve()
    medoids_p = Path(args.medoids_json).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not struct_p.exists() or not medoids_p.exists():
        raise FileNotFoundError("FATAL: Required reference structure or medoids summary missing!")

    clean_xtcs = []
    for raw_arg in args.xtc:
        for chunk in str(raw_arg).split():
            c_clean = chunk.strip()
            if c_clean:
                p_xtc = Path(c_clean).resolve()
                if p_xtc.exists() and p_xtc.stat().st_size > 0:
                    clean_xtcs.append(str(p_xtc))

    if not clean_xtcs:
        raise ValueError("FATAL: No valid trajectory files provided to extract_centroids!")

    with open(medoids_p, "r", encoding="utf-8") as f:
        medoids = json.load(f)

    u = mda.Universe(str(struct_p), clean_xtcs)
    protein_atoms = u.select_atoms("protein")
    if len(protein_atoms) == 0:
        raise ValueError("FATAL: No protein atoms found in reference structure!")

    # Discover all small molecules (active ligand, cofactors, ions) to preserve system integrity
    small_mol_resnames = set()
    for res in u.residues:
        r_name = res.resname.strip().upper()
        if r_name not in STANDARD_PROTEIN_RESIDUES and r_name not in WATER_RESIDUES:
            small_mol_resnames.add(r_name)

    active_res = args.ligand_resname.strip().upper() if args.ligand_resname else ""

    # Parse heterogens specification to obtain templates for ALL organic molecules (ligand and cofactors)
    try:
        hetero_specs = json.loads(args.heterogens_spec) if args.heterogens_spec else []
    except Exception as e:
        raise ValueError(f"FATAL: Failed to parse heterogens specification JSON: {e}")

    hetero_by_res = {item["res_name"].strip().upper(): item for item in hetero_specs if "res_name" in item}
    tmpl_dir_p = Path(args.templates_dir).resolve() if args.templates_dir else None

    # Pre-load canonical template molecules for all declared organic small molecules
    templates_by_res = {}
    for r_code, spec in hetero_by_res.items():
        if spec.get("is_ion", False):
            continue  # Monoatomic ions have no bond-order templates

        orig_sdf_str = spec.get("template_sdf", "")
        cand_paths = []
        if orig_sdf_str:
            cand_paths.append(Path(orig_sdf_str).resolve())
            cand_paths.append(Path(Path(orig_sdf_str).name).resolve())
            if tmpl_dir_p and tmpl_dir_p.exists():
                cand_paths.append(tmpl_dir_p / Path(orig_sdf_str).name)

        found_mol = None
        for p in cand_paths:
            if p.exists() and p.stat().st_size > 0:
                suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=False)
                found_mol = next((m for m in suppl if m is not None), None)
                if found_mol:
                    break

        if found_mol is not None:
            templates_by_res[r_code] = found_mol
            print(f"[EXTRACT_CENTROIDS] Loaded pristine topology template for small molecule '{r_code}'.", flush=True)

    # Sort medoid frames monotonically for optimal linear seeking
    sorted_medoids = sorted(medoids, key=lambda x: x["medoid_frame"])
    manifest_by_rank = {}
    batch_convert_tasks = []

    for stat in sorted_medoids:
        rank = stat["rank"]
        m_frame = stat["medoid_frame"]
        pop_pct = stat["population_pct"]

        # Monotonic seek forward in trajectory
        u.trajectory[m_frame]

        # 1. Export pure apo protein receptor
        rec_pdb = outdir / f"centroid_{rank}_receptor.pdb"
        rec_mae = outdir / f"centroid_{rank}_receptor.mae"
        protein_atoms.write(str(rec_pdb))
        batch_convert_tasks.append({"in_path": str(rec_pdb), "out_path": str(rec_mae)})

        # 2. Export full assembled complex view
        view_pdb = outdir / f"centroid_{rank}_complex.pdb"
        view_mae = outdir / f"centroid_{rank}_complex.mae"
        u.atoms.write(str(view_pdb))
        batch_convert_tasks.append({"in_path": str(view_pdb), "out_path": str(view_mae)})

        extracted_poses = {}
        # Iterate over individual small molecule residues to support multiple copies (ions, oligomeric ligands)
        for res in u.residues:
            r_name = res.resname.strip().upper()
            if r_name not in small_mol_resnames:
                continue

            mol_atoms = res.atoms
            if len(mol_atoms) == 0:
                continue

            # Deterministic residue suffix if multiple copies of the same small molecule exist in the system
            same_res_count = len(u.select_atoms(f"resname {r_name}").residues)
            res_key = f"{r_name}_{res.resid}" if same_res_count > 1 else r_name
            dest_sdf = outdir / f"centroid_{rank}_{res_key}_pose.sdf"
            dest_mae = outdir / f"centroid_{rank}_{res_key}_pose.mae"

            # Strict in-memory PDB serialization preserving buffer on writer close
            pdb_stream = NonClosingStringIO()
            with mda.coordinates.PDB.PDBWriter(pdb_stream) as writer:
                writer.write(mol_atoms)
            raw_pdb_block = pdb_stream.getvalue()

            raw_mol = Chem.MolFromPDBBlock(raw_pdb_block, removeHs=False, sanitize=False)
            if not raw_mol:
                raise ValueError(f"FATAL: Failed to parse small molecule '{res_key}' from MD coordinates at frame {m_frame}!")

            spec = hetero_by_res.get(r_name)
            is_active_ligand = spec.get("is_ligand", (r_name == active_res)) if spec else (r_name == active_res)
            is_ion = spec.get("is_ion", False) if spec else False
            tmpl = templates_by_res.get(r_name)

            # Rigorous topological transfer with full automorphism symmetry check
            # APPLIES IDENTICALLY TO LIGAND AND ALL ORGANIC COFACTORS!
            if tmpl is not None and not is_ion:
                restored_pose = transfer_coordinates_to_template(
                    template_mol=tmpl,
                    source_mol=raw_mol,
                    relax_missing_hydrogens=False,
                )
            else:
                # Monoatomic ions (MG, ZN, CA, CL) retain raw Cartesian coordinates
                restored_pose = raw_mol

            with Chem.SDWriter(str(dest_sdf)) as w:
                w.SetForceV3000(False)
                w.write(restored_pose)

            batch_convert_tasks.append({"in_path": str(dest_sdf), "out_path": str(dest_mae)})
            extracted_poses[res_key] = {
                "res_name": r_name,
                "resid": int(res.resid),
                "is_ligand": is_active_ligand,
                "is_ion": is_ion,
                "sdf_path": dest_sdf.name,
                "mae_path": dest_mae.name,
            }

        manifest_by_rank[rank] = {
            "rank": rank,
            "population_pct": pop_pct,
            "medoid_frame": m_frame,
            "medoid_time_ns": stat["medoid_time_ns"],
            "medoid_replicate": stat["medoid_replicate"],
            "receptor_pdb": rec_pdb.name,
            "receptor_mae": rec_mae.name,
            "complex_pdb": view_pdb.name,
            "complex_mae": view_mae.name,
            "poses": extracted_poses,
        }
        print(f"[EXTRACT_CENTROIDS] Centroid #{rank} extracted (frame {m_frame}, {pop_pct}%).", flush=True)

    # Restore canonical rank sorting in final manifest
    final_manifest = [manifest_by_rank[r] for r in sorted(manifest_by_rank.keys())]

    with open(outdir / "extracted_manifest.json", "w", encoding="utf-8") as f:
        json.dump(final_manifest, f, indent=2)

    with open(outdir / "batch_convert_tasks.json", "w", encoding="utf-8") as f:
        json.dump(batch_convert_tasks, f, indent=2)

    print(f"✓ Extracted {len(final_manifest)} centroid structures and emitted batch conversion manifest.", flush=True)


if __name__ == "__main__":
    main()
