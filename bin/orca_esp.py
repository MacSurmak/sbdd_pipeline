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
    parser = argparse.ArgumentParser(
        description="ORCA Single-Point ESP Driver via Archetype Template"
    )
    parser.add_argument(
        "--template", required=True, help="Path to ORCA .inp template file"
    )
    parser.add_argument("--sdf", required=True, help="Input relaxed molecule structure")
    parser.add_argument("--charge", required=True, type=int, help="Net formal charge")
    parser.add_argument("--mult", required=True, type=int, help="Spin multiplicity")
    parser.add_argument("--nprocs", type=int, required=True, help="Number of CPU cores")
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

    tmpl_path = Path(args.template).resolve()
    if not tmpl_path.exists() or tmpl_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: ORCA template file does not exist or is empty: {tmpl_path}"
        )

    suppl = Chem.SDMolSupplier(str(in_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(
            f"FATAL: Could not read valid molecular structure from {in_path}"
        )

    conf = mol.GetConformer()
    n_atoms = mol.GetNumAtoms()

    coord_block = "\n".join(
        f"  {mol.GetAtomWithIdx(idx).GetSymbol():<3s} {conf.GetAtomPosition(idx).x:12.6f} {conf.GetAtomPosition(idx).y:12.6f} {conf.GetAtomPosition(idx).z:12.6f}"
        for idx in range(n_atoms)
    )
    tmpl_text = tmpl_path.read_text(encoding="utf-8")
    inp_content = (
        tmpl_text.replace("{NPROCS}", str(args.nprocs))
        .replace("{MAXCORE}", str(args.maxcore_mb))
        .replace("{CHARGE}", str(args.charge))
        .replace("{MULT}", str(args.mult))
        .replace("{COORDINATES}", coord_block)
    )

    inp_file = Path("orca_esp.inp")
    inp_file.write_text(inp_content, encoding="utf-8")
    out_file = Path(args.out_log).resolve()

    run_orca(inp_file, out_file, args.nprocs)
    print(
        f"✓ Successfully calculated ORCA ESP grid ({tmpl_path.stem}): {out_file.name}",
        flush=True,
    )


if __name__ == "__main__":
    main()
