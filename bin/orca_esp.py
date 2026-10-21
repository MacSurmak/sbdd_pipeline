#!/usr/bin/env python3
"""
bin/orca_esp.py
Atomic ORCA Electrostatic Potential (ESP) single-point calculation driver.
Evaluates canonical Merz-Kollman grid in vacuum with TightSCF.
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
    """Executes ORCA single-point calculation with real-time stdout streaming."""
    orca_bin = shutil.which("orca") or "/opt/containers/bin/orca"
    env = setup_orca_environment(nprocs)

    print(
        f"Exec: {orca_bin} {inp_file.name} (live streaming to .command.log)",
        flush=True,
    )
    with open(out_file, "w", encoding="utf-8") as out_f:
        proc = subprocess.Popen(
            [orca_bin, inp_file.name],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
        )
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            out_f.write(line)
            out_f.flush()
        proc.wait()

    if not out_file.exists() or out_file.stat().st_size == 0:
        raise RuntimeError(f"FATAL: ORCA failed to produce output log: {out_file}")

    log_content = out_file.read_text(encoding="utf-8", errors="ignore")
    if proc.returncode != 0 or "ORCA TERMINATED NORMALLY" not in log_content:
        tail_lines = "\n".join(log_content.splitlines()[-40:])
        raise RuntimeError(
            f"FATAL: ORCA ESP calculation terminated abnormally (exit code {proc.returncode})!\n"
            f"--- TAIL OF {out_file.name} ---\n{tail_lines}"
        )


def main():
    parser = argparse.ArgumentParser(description="ORCA Single-Point ESP Driver")
    parser.add_argument(
        "--sdf", required=True, help="Input relaxed molecule structure"
    )
    parser.add_argument("--charge", required=True, type=int, help="Net formal charge")
    parser.add_argument("--mult", required=True, type=int, help="Spin multiplicity")
    parser.add_argument(
        "--theory", default="HF", help="ESP Theory level (canonical: HF)"
    )
    parser.add_argument(
        "--basis", default="6-31G*", help="ESP Basis set (canonical: 6-31G*)"
    )
    parser.add_argument(
        "--density", default="Merz-Kollman", help="ESP density algorithm"
    )
    parser.add_argument(
        "--scf_convergence", default="TightSCF", help="SCF convergence criterion"
    )
    parser.add_argument(
        "--scf_maxiter", type=int, default=250, help="Max SCF iterations"
    )
    parser.add_argument(
        "--nprocs", type=int, required=True, help="Number of CPU cores"
    )
    parser.add_argument(
        "--maxcore_mb", type=int, required=True, help="Per-core RAM allocation in MB"
    )
    parser.add_argument(
        "--out_log", default="orca_esp.out", help="Output execution log"
    )
    args = parser.parse_args()

    in_path = Path(args.sdf).resolve()
    if not in_path.exists() or in_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Input SDF deliverable not found: {in_path}")

    suppl = Chem.SDMolSupplier(str(in_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(
            f"FATAL: Could not read valid molecular structure from {in_path}"
        )

    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()

    inp_lines = [
        f"! {args.theory} {args.basis} RESP {args.scf_convergence}",
        f"%maxcore {args.maxcore_mb}",
        f"%pal nprocs {args.nprocs} end",
        f"%scf maxiter {args.scf_maxiter} end",
        f"* xyz {args.charge} {args.mult}",
    ]
    for idx in range(n_atoms):
        pos = conf.GetAtomPosition(idx)
        elem = mol.GetAtomWithIdx(idx).GetSymbol()
        inp_lines.append(f"  {elem:<3s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}")
    inp_lines.append("*\n")

    inp_file = Path("orca_esp.inp")
    inp_file.write_text("\n".join(inp_lines), encoding="utf-8")
    out_file = Path(args.out_log).resolve()

    run_orca(inp_file, out_file, args.nprocs)
    print(
        f"✓ Successfully calculated ORCA ESP grid ({args.theory}/{args.basis}): {out_file.name}",
        flush=True,
    )


if __name__ == "__main__":
    main()
