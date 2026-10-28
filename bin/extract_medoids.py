#!/usr/bin/env python3
"""
bin/extract_medoids.py
Slurm Worker Script: Extracts physical medoid complexes (Protein, Ligand, Cofactors/Ions)
from MD trajectories in $SLURM_TMPDIR without PDB element corruption.
Produces PDB v3.3, converts to MAE via structconvert, and emits structured TSV receipt.
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
import MDAnalysis as mda
from rdkit import Chem

# Inject framework root for canonical chem engine
FRAMEWORK_DIR = Path(__file__).resolve().parent.parent
if str(FRAMEWORK_DIR) not in sys.path:
    sys.path.insert(0, str(FRAMEWORK_DIR))

from cli.core.chem import transfer_coordinates_to_template


def convert_to_mae(src_path: Path, dst_path: Path):
    """Converts PDB/SDF to MAE using Schrödinger structconvert with strict Fail-Fast enforcement."""
    structconvert_bin = shutil.which("structconvert") or "/opt/containers/bin/structconvert"
    if Path(structconvert_bin).exists() or shutil.which("structconvert"):
        cmd = [structconvert_bin, str(src_path), str(dst_path)]
    elif shutil.which("schrodinger"):
        cmd = ["schrodinger", "run", "structconvert", str(src_path), str(dst_path)]
    else:
        raise RuntimeError(
            "FATAL: Neither 'structconvert' nor 'schrodinger' binary found in PATH!"
        )

    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0 or not dst_path.exists() or dst_path.stat().st_size == 0:
        raise RuntimeError(
            f"FATAL: structconvert failed for {src_path} -> {dst_path} (exit {res.returncode}):\n"
            f"{res.stderr}\n{res.stdout}"
        )


def main():
    parser = argparse.ArgumentParser(description="Extract medoid complexes on Slurm node")
    parser.add_argument("--spec", required=True, help="Task specification JSON")
    parser.add_argument("--outdir", default=".", help="Local output scratch directory")
    args = parser.parse_args()

    spec_file = Path(args.spec).resolve()
    if not spec_file.exists():
        raise FileNotFoundError(f"FATAL: Specification JSON not found: {spec_file}")

    with open(spec_file, "r", encoding="utf-8") as f:
        spec = json.load(f)

    target_name = spec["target_name"]
    system_name = spec["system_name"]
    clustering_name = spec["clustering_name"]
    parent_frame_name = spec["parent_frame_name"]
    parent_frame_id = spec["parent_frame_id"]
    feature_group = spec["feature_group"]
    align_group = spec["align_group"]
    dim_reduction = spec["dim_reduction"]

    ref_pdb_path = Path(spec["ref_pdb"]).resolve()
    clean_xtcs = [str(Path(x).resolve()) for x in spec["clean_xtcs"]]

    if not ref_pdb_path.exists():
        raise FileNotFoundError(f"FATAL: Reference PDB missing: {ref_pdb_path}")

    # Load trajectory and strictly guess missing topology elements to guarantee PDB v3.3 standard compliance
    u = mda.Universe(str(ref_pdb_path), clean_xtcs)
    u.guess_TopologyAttrs(to_guess=["elements", "masses"])

    protein_atoms = u.select_atoms("protein")
    if len(protein_atoms) == 0:
        raise ValueError(f"FATAL: Zero protein atoms found in {ref_pdb_path.name}!")

    lig_spec = spec.get("ligand")
    cofactors_spec = spec.get("cofactors", [])
    clusters = spec.get("clusters", [])

    out_base = Path(args.outdir).resolve()
    receipt_lines = []

    for cl in clusters:
        rank = cl["rank"]
        m_frame = cl["medoid_frame"]
        pop = cl["population"]
        spread = cl["spread"]

        # Monotonic trajectory seek
        u.trajectory[m_frame]

        new_frame_name = f"{parent_frame_name}_{clustering_name}_c{rank}"
        frame_dir = out_base / "targets" / target_name / "frames" / new_frame_name
        frame_dir.mkdir(parents=True, exist_ok=True)

        rec_pdb = frame_dir / "receptor.pdb"
        rec_mae = frame_dir / "receptor.mae"
        protein_atoms.write(str(rec_pdb))
        convert_to_mae(rec_pdb, rec_mae)

        rel_rec_pdb = f"targets/{target_name}/frames/{new_frame_name}/receptor.pdb"
        rel_rec_mae = f"targets/{target_name}/frames/{new_frame_name}/receptor.mae" if rec_mae.exists() and rec_mae.stat().st_size > 0 else "None"

        receipt_lines.append(
            f"FRAME\t{target_name}\t{new_frame_name}\t{rel_rec_pdb}\t{rel_rec_mae}\t"
            f"{parent_frame_id}\t{rank}\t{pop}\t{spread}\t{feature_group}\t{align_group}\t{dim_reduction}"
        )

        def _extract_small_molecule(mol_spec: dict, is_primary_ligand: bool) -> tuple[str, str]:
            res_clean = mol_spec["res_name"].strip().upper()
            inchikey = mol_spec["inchikey"]
            state_id = mol_spec["state_id"]
            is_ion = mol_spec.get("is_ion", False)

            mol_atoms = u.select_atoms(f"resname {res_clean}")
            if len(mol_atoms) == 0:
                raise ValueError(f"FATAL: Molecule '{res_clean}' not found in frame {m_frame} of {system_name}!")

            tmpl_path = Path(mol_spec["template_sdf"]).resolve()
            suppl = Chem.SDMolSupplier(str(tmpl_path), removeHs=False, sanitize=False)
            tmpl_mol = next((m for m in suppl if m is not None), None)
            if tmpl_mol is None:
                raise ValueError(f"FATAL: Failed to read reference SDF template from {tmpl_path}")

            raw_positions = mol_atoms.positions
            if is_ion:
                restored_pose = Chem.Mol(tmpl_mol)
                restored_pose.GetConformer().SetAtomPosition(
                    0, [float(raw_positions[0][0]), float(raw_positions[0][1]), float(raw_positions[0][2])]
                )
            else:
                mda_elements = []
                for a in mol_atoms:
                    if hasattr(a, "element") and a.element and a.element.strip():
                        mda_elements.append(a.element.strip().capitalize())
                    else:
                        clean_n = "".join([c for c in a.name if not c.isdigit()]).strip()
                        mda_elements.append(
                            clean_n[:2].capitalize()
                            if len(clean_n) >= 2 and clean_n[:2].capitalize() in {"Cl", "Br", "Fe", "Mg", "Zn", "Ca", "Na"}
                            else clean_n[:1].upper()
                        )

                restored_pose = transfer_coordinates_to_template(
                    template_mol=tmpl_mol,
                    source_mol_or_coords=raw_positions,
                    source_elements=mda_elements,
                    relax_missing_hydrogens=False,
                )

            pose_name = f"{new_frame_name}_{res_clean}_centroid"
            pose_dir = out_base / "substances" / inchikey / "poses" / pose_name
            pose_dir.mkdir(parents=True, exist_ok=True)
            pose_sdf = pose_dir / f"{res_clean}_pose.sdf"
            pose_mae = pose_dir / f"{res_clean}_pose.mae"

            with Chem.SDWriter(str(pose_sdf)) as w:
                w.SetForceV3000(False)
                w.write(restored_pose)
            convert_to_mae(pose_sdf, pose_mae)

            rel_sdf = f"substances/{inchikey}/poses/{pose_name}/{res_clean}_pose.sdf"
            rel_mae = f"substances/{inchikey}/poses/{pose_name}/{res_clean}_pose.mae" if pose_mae.exists() and pose_mae.stat().st_size > 0 else "None"

            lig_flag = "1" if is_primary_ligand else "0"
            receipt_lines.append(
                f"POSE\t{new_frame_name}\t{pose_name}\t{inchikey}\t{res_clean}\t{state_id}\t{rel_sdf}\t{rel_mae}\t{lig_flag}"
            )
            return pose_name, res_clean

        lig_pose_name = ""
        lig_res_name = ""
        if lig_spec:
            lig_pose_name, lig_res_name = _extract_small_molecule(lig_spec, is_primary_ligand=True)

        cof_pose_names = []
        for cof in cofactors_spec:
            c_pname, _ = _extract_small_molecule(cof, is_primary_ligand=False)
            cof_pose_names.append(c_pname)

        sys_name = f"{new_frame_name}_{lig_res_name}" if lig_res_name else new_frame_name
        cof_str = ",".join(cof_pose_names)
        receipt_lines.append(f"SYSTEM\t{target_name}\t{sys_name}\t{new_frame_name}\t{lig_pose_name}\t{cof_str}")

    out_receipt = out_base / "extract_medoids_receipt.tsv"
    out_receipt.write_text("\n".join(receipt_lines) + "\n", encoding="utf-8")
    print(f"[EXTRACT_MEDOIDS] Successfully extracted {len(clusters)} complexes for {system_name}.", flush=True)


if __name__ == "__main__":
    main()
