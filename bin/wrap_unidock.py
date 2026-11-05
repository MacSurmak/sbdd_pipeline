#!/usr/bin/env python3
"""
bin/wrap_unidock.py
High-throughput GPU docking worker for Uni-Dock.
Executes entirely in node-local scratch ($SLURM_TMPDIR), computes pocket box via SSOT
cli.core.chem.compute_autobox_geometry, runs batch GPU sampling, and restores exact
chemical topology and explicit hydrogens via transfer_coordinates_to_template.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

FRAMEWORK_DIR = Path(__file__).resolve().parent.parent
if str(FRAMEWORK_DIR) not in sys.path:
    sys.path.insert(0, str(FRAMEWORK_DIR))

from cli.core.chem import (
    compute_autobox_geometry,
    convert_receptor_pdb_to_pdbqt,
    transfer_coordinates_to_template,
)


def _extract_unidock_score(mol: Chem.Mol) -> float:
    """Extracts Vina affinity (kcal/mol) from Uni-Dock output SDF properties."""
    for prop_name in mol.GetPropNames():
        p_clean = prop_name.strip()
        if p_clean.lower() in ("uni-dock result", "energy", "affinity", "minimizedaffinity", "vina result"):
            raw = mol.GetProp(prop_name)
            m = re.search(r"ENERGY\s*=\s*([-+]?\d+(?:\.\d+)?)", raw, re.IGNORECASE)
            if m:
                return float(m.group(1))
            m = re.search(r"[-+]?\d+(?:\.\d+)?", raw)
            if m:
                return float(m.group(0))
    raise ValueError(
        f"FATAL: Uni-Dock output pose '{mol.GetProp('_Name') if mol.HasProp('_Name') else '?'}' "
        f"is missing energy property! Available properties: {list(mol.GetPropNames())}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Uni-Dock GPU batch docking worker")
    parser.add_argument("--receptor", required=True, help="Input receptor (.pdbqt or .pdb)")
    parser.add_argument("--ligand", required=True, help="Input batch ligands (.sdf)")
    parser.add_argument("--ref_ligand", required=True, help="Reference pocket ligand (.sdf)")
    parser.add_argument("--search_mode", required=True, choices=["fast", "balance", "detail"])
    parser.add_argument("--scoring", required=True, choices=["vina", "vinardo", "ad4"])
    parser.add_argument("--refine_step", type=int, required=True)
    parser.add_argument("--num_modes", type=int, required=True)
    parser.add_argument("--energy_range", type=float, required=True)
    parser.add_argument("--autobox_add", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out_sdf", required=True, help="Output merged docked poses (.sdf)")
    parser.add_argument("--out_log", required=True, help="Output execution log")
    args = parser.parse_args()

    rec_path = Path(args.receptor).resolve()
    lig_path = Path(args.ligand).resolve()
    ref_path = Path(args.ref_ligand).resolve()

    if not rec_path.exists() or rec_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Receptor file missing or empty: {rec_path}")
    if not lig_path.exists() or lig_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Ligand chunk SDF missing or empty: {lig_path}")

    box_geom = compute_autobox_geometry(ref_path, args.autobox_add)
    cx, cy, cz = box_geom["autobox_center"]
    sx, sy, sz = box_geom["autobox_size"]

    scratch_root = os.environ.get("SLURM_TMPDIR") or "/tmp"
    with tempfile.TemporaryDirectory(prefix="unidock_work_", dir=scratch_root) as work_str:
        work_dir = Path(work_str)
        in_dir = work_dir / "in_ligs"
        out_dir = work_dir / "out_ligs"
        in_dir.mkdir(parents=True, exist_ok=True)
        out_dir.mkdir(parents=True, exist_ok=True)

        if rec_path.suffix.lower() == ".pdbqt":
            rec_pdbqt = rec_path
        else:
            rec_pdbqt = work_dir / "receptor.pdbqt"
            convert_receptor_pdb_to_pdbqt(rec_path, rec_pdbqt)

        templates_by_idx: dict[int, tuple[str, Chem.Mol]] = {}
        index_paths: list[str] = []

        for idx, tmpl_mol in enumerate(
            Chem.SDMolSupplier(str(lig_path), removeHs=False, sanitize=False)
        ):
            if tmpl_mol is None:
                raise ValueError(
                    f"FATAL: Corrupted molecule record at index {idx} in {lig_path}!"
                )
            t_name = (
                tmpl_mol.GetProp("_Name").strip()
                if tmpl_mol.HasProp("_Name")
                else ""
            )
            if not t_name and tmpl_mol.HasProp("s_m_title"):
                t_name = tmpl_mol.GetProp("s_m_title").strip()
            if not t_name:
                raise ValueError(
                    f"FATAL: Unnamed molecule at index {idx} in {lig_path}!"
                )
            tmpl_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(tmpl_mol)
            templates_by_idx[idx] = (t_name, Chem.Mol(tmpl_mol))

            single_sdf = in_dir / f"lig_{idx:05d}.sdf"
            mol_copy = Chem.Mol(tmpl_mol)
            mol_copy.SetProp("_Name", f"{t_name}__tmpl_{idx}")
            with Chem.SDWriter(str(single_sdf)) as w:
                w.SetForceV3000(False)
                w.write(mol_copy)
            index_paths.append(str(single_sdf))

        if not index_paths:
            raise ValueError(f"FATAL: Zero valid molecules found in {lig_path}!")

        ligand_index_file = work_dir / "ligand_index.txt"
        ligand_index_file.write_text("\n".join(index_paths) + "\n", encoding="utf-8")

        unidock_bin = shutil.which("unidock") or "/opt/containers/bin/unidock"
        cmd = [
            unidock_bin,
            "--receptor",
            str(rec_pdbqt),
            "--ligand_index",
            str(ligand_index_file),
            "--center_x",
            f"{cx:.4f}",
            "--center_y",
            f"{cy:.4f}",
            "--center_z",
            f"{cz:.4f}",
            "--size_x",
            f"{sx:.4f}",
            "--size_y",
            f"{sy:.4f}",
            "--size_z",
            f"{sz:.4f}",
            "--search_mode",
            args.search_mode,
            "--scoring",
            args.scoring,
            "--refine_step",
            str(args.refine_step),
            "--num_modes",
            str(args.num_modes),
            "--energy_range",
            str(args.energy_range),
            "--seed",
            str(args.seed),
            "--dir",
            str(out_dir),
        ]

        print(f"Exec: {' '.join(cmd)}", flush=True)
        with open(args.out_log, "w", encoding="utf-8") as log_f:
            res = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, check=False)

        if res.returncode != 0:
            log_text = Path(args.out_log).read_text(errors="ignore") if Path(args.out_log).exists() else ""
            print(
                f"FATAL: Uni-Dock failed with exit code {res.returncode}!\n"
                f"Command: {' '.join(cmd)}\n"
                f"==================== UNIDOCK LOG ({args.out_log}) ====================\n"
                f"{log_text}\n"
                f"======================================================================",
                file=sys.stderr,
                flush=True,
            )
            sys.exit(res.returncode)

        # Group restored poses by substance title and sort modes strictly by Vina affinity
        poses_by_title: dict[str, list[tuple[float, Chem.Mol]]] = {}
        for idx in sorted(templates_by_idx.keys()):
            orig_title, exact_tmpl = templates_by_idx[idx]
            out_lig_sdf = out_dir / f"lig_{idx:05d}_out.sdf"
            if not out_lig_sdf.exists() or out_lig_sdf.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Uni-Dock failed to produce output SDF for record {idx} ('{orig_title}'): {out_lig_sdf}"
                )

            suppl = Chem.SDMolSupplier(str(out_lig_sdf), removeHs=False, sanitize=False)
            for pose_mol in suppl:
                if pose_mol is None:
                    continue
                vina_score = _extract_unidock_score(pose_mol)
                restored_mol = transfer_coordinates_to_template(
                    template_mol=exact_tmpl,
                    source_mol=pose_mol,
                    relax_missing_hydrogens=False,
                )
                restored_mol.SetProp("_Name", orig_title)
                restored_mol.SetProp("s_m_title", orig_title)
                restored_mol.SetProp("unidock_vina_score", f"{vina_score:.3f}")
                restored_mol.SetProp("minimizedAffinity", f"{vina_score:.3f}")
                poses_by_title.setdefault(orig_title, []).append((vina_score, restored_mol))

        tmp_out_sdf = f"{args.out_sdf}.tmp"
        total_written = 0
        with Chem.SDWriter(tmp_out_sdf) as writer:
            writer.SetForceV3000(False)
            for orig_title, pose_list in poses_by_title.items():
                # Sort modes across all protomers of the substance by best (lowest) Vina energy
                pose_list.sort(key=lambda x: x[0])
                for rank_idx, (_, r_mol) in enumerate(pose_list[: args.num_modes], start=1):
                    r_mol.SetProp("unidock_mode_rank", str(rank_idx))
                    writer.write(r_mol)
                    total_written += 1

        if total_written == 0:
            raise RuntimeError(
                f"FATAL: Zero valid docked poses written by wrap_unidock.py for {lig_path}!"
            )

        os.replace(tmp_out_sdf, args.out_sdf)
        print(
            f"✓ Uni-Dock completed: {total_written} topology-validated poses across {len(poses_by_title)} substance(s) -> {args.out_sdf}",
            flush=True,
        )


if __name__ == "__main__":
    main()
