#!/usr/bin/env python3
"""
bin/orca_opt.py
Atomic ORCA geometry optimization driver.
Executes DFT relaxation (pure GGA pre-opt or hybrid refinement) with implicit solvent.
Strict Fail-Fast on unconverged SCF or abnormal termination.
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from rdkit import Chem


def setup_orca_environment(nprocs: int) -> dict:
    """Prepares execution environment for OpenMPI/PRRTE under Slurm cgroups."""
    env = os.environ.copy()
    cpus_per_task = env.get("SLURM_CPUS_PER_TASK", str(nprocs))
    env["SLURM_TASKS_PER_NODE"] = str(cpus_per_task)
    env["SLURM_NTASKS"] = str(cpus_per_task)
    env["SLURM_NPROCS"] = str(cpus_per_task)

    env["PRTE_MCA_rmaps_default_mapping_policy"] = ":oversubscribe"
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    env["OMPI_MCA_hwloc_base_binding_policy"] = "none"
    env["OMPI_MCA_btl"] = "^openib"
    return env


def run_orca(inp_file: Path, out_file: Path, nprocs: int):
    """Executes ORCA binary with strict error trapping."""
    orca_bin = shutil.which("orca") or "/opt/containers/bin/orca"
    env = setup_orca_environment(nprocs)

    print(f"Exec: {orca_bin} {inp_file.name} > {out_file.name}", flush=True)
    with open(out_file, "w", encoding="utf-8") as out_f:
        res = subprocess.run(
            [orca_bin, inp_file.name],
            stdout=out_f,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            env=env,
        )

    if not out_file.exists() or out_file.stat().st_size == 0:
        raise RuntimeError(f"FATAL: ORCA failed to produce output log: {out_file}")

    log_content = out_file.read_text(encoding="utf-8", errors="ignore")
    if res.returncode != 0 or "ORCA TERMINATED NORMALLY" not in log_content:
        tail_lines = "\n".join(log_content.splitlines()[-40:])
        raise RuntimeError(
            f"FATAL: ORCA terminated abnormally with exit code {res.returncode}!\n"
            f"--- TAIL OF {out_file.name} ---\n{tail_lines}"
        )


def main():
    parser = argparse.ArgumentParser(description="ORCA Atomic Geometry Optimizer")
    parser.add_argument("--input", required=True, help="Input molecule (.sdf or .xyz)")
    parser.add_argument("--charge", required=True, type=int, help="Formal charge")
    parser.add_argument("--mult", required=True, type=int, help="Spin multiplicity")
    parser.add_argument("--theory", required=True, help="Functional (e.g. PBE, B3LYP)")
    parser.add_argument("--basis", required=True, help="Basis set (e.g. def2-SVP)")
    parser.add_argument("--dispersion", default="", help="Dispersion (e.g. D4)")
    parser.add_argument("--solvent", default="", help="Solvent (e.g. CPCM(Water))")
    parser.add_argument(
        "--rijcosx", action="store_true", help="Enable RIJCOSX approximation"
    )
    parser.add_argument(
        "--maxiter", type=int, required=True, help="Max geometry cycles"
    )
    parser.add_argument(
        "--scf_convergence", required=True, help="SCF convergence (e.g. NormalSCF)"
    )
    parser.add_argument(
        "--scf_maxiter", type=int, required=True, help="Max SCF cycles"
    )
    parser.add_argument(
        "--nprocs", type=int, required=True, help="Processors (%pal nprocs)"
    )
    parser.add_argument(
        "--maxcore_mb", type=int, required=True, help="Per-core RAM in MB"
    )
    parser.add_argument("--out_xyz", default="", help="Output optimized XYZ file")
    parser.add_argument("--out_sdf", default="", help="Output optimized SDF file")
    parser.add_argument(
        "--out_log", default="orca_opt.out", help="Output execution log"
    )
    args = parser.parse_args()

    in_path = Path(args.input).resolve()
    if not in_path.exists():
        raise FileNotFoundError(f"FATAL: Input file does not exist: {in_path}")

    atoms_data = []
    mol_rdkit = None

    if in_path.suffix.lower() == ".sdf":
        suppl = Chem.SDMolSupplier(str(in_path), removeHs=False, sanitize=False)
        mol_rdkit = next((m for m in suppl if m is not None), None)
        if mol_rdkit is None:
            raise ValueError(f"FATAL: Failed to parse RDKit molecule from {in_path}")
        conf = mol_rdkit.GetConformer()
        for idx in range(mol_rdkit.GetNumAtoms()):
            pos = conf.GetAtomPosition(idx)
            elem = mol_rdkit.GetAtomWithIdx(idx).GetSymbol()
            atoms_data.append((elem, pos.x, pos.y, pos.z))
    elif in_path.suffix.lower() == ".xyz":
        lines = in_path.read_text(encoding="utf-8", errors="ignore").splitlines()
        coord_lines = [l.strip() for l in lines if l.strip() and not l.startswith("#")][
            2:
        ]
        for line in coord_lines:
            parts = line.split()
            atoms_data.append(
                (parts[0], float(parts[1]), float(parts[2]), float(parts[3]))
            )
    else:
        raise ValueError(f"FATAL: Unsupported input format for ORCA: {in_path.suffix}")

    # Build calculation header
    rij_flag = "RIJCOSX def2/J" if args.rijcosx else "RI def2/J"
    disp_flag = args.dispersion if args.dispersion.lower() != "none" else ""
    solv_flag = args.solvent if args.solvent.lower() != "none" else ""

    inp_lines = [
        f"! {args.theory} {disp_flag} {args.basis} {rij_flag} {solv_flag} Opt {args.scf_convergence}",
        f"%maxcore {args.maxcore_mb}",
        f"%pal nprocs {args.nprocs} end",
        f"%geom maxiter {args.maxiter} end",
        f"%scf maxiter {args.scf_maxiter} end",
        f"* xyz {args.charge} {args.mult}",
    ]
    for elem, x, y, z in atoms_data:
        inp_lines.append(f"  {elem:<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
    inp_lines.append("*\n")

    inp_file = Path("orca_opt.inp")
    inp_file.write_text("\n".join(inp_lines), encoding="utf-8")
    out_file = Path(args.out_log)

    run_orca(inp_file, out_file, args.nprocs)

    orca_xyz = Path("orca_opt.xyz")
    if not orca_xyz.exists() or orca_xyz.stat().st_size == 0:
        raise FileNotFoundError(
            "FATAL: ORCA completed but 'orca_opt.xyz' was not generated! "
            "Optimization did not converge."
        )

    # Read optimized coordinates
    opt_lines = orca_xyz.read_text(encoding="utf-8", errors="ignore").splitlines()
    opt_coords = [l.strip() for l in opt_lines if l.strip() and not l.startswith("#")][
        2:
    ]

    if len(opt_coords) != len(atoms_data):
        raise ValueError(
            f"FATAL: Optimized atom count ({len(opt_coords)}) does not match input ({len(atoms_data)})!"
        )

    if args.out_xyz:
        shutil.copy2(orca_xyz, Path(args.out_xyz).resolve())

    if args.out_sdf:
        if mol_rdkit is None:
            raise ValueError(
                "FATAL: Cannot export SDF deliverable when input was pure XYZ."
            )
        conf = mol_rdkit.GetConformer()
        for idx, line in enumerate(opt_coords):
            parts = line.split()
            conf.SetAtomPosition(
                idx, (float(parts[1]), float(parts[2]), float(parts[3]))
            )
        out_sdf_path = Path(args.out_sdf).resolve()
        out_sdf_path.parent.mkdir(parents=True, exist_ok=True)
        with Chem.SDWriter(str(out_sdf_path)) as writer:
            writer.SetForceV3000(False)
            writer.write(mol_rdkit)

    print(
        f"✓ Successfully converged ORCA optimization ({args.theory}): {len(opt_coords)} atoms.",
        flush=True,
    )


if __name__ == "__main__":
    main()
