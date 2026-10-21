#!/usr/bin/env python3
"""
bin/build_gaff2_topology.py
Generates GAFF2 mol2 and frcmod parameter files using AmberTools antechamber and parmchk2.
Preserves the pristine coordinates of the input reference SDF while applying external QM charges (-c rc -cf).
Strict Fail-Fast error trapping.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from rdkit import Chem


def run_cmd(cmd: list[str], desc: str = "Command"):
    """Executes an external command with strict Fail-Fast verification."""
    print(f"[{desc}] Exec: {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        print(
            f"FATAL: {desc} failed with exit code {res.returncode}!\n"
            f"--- STDOUT ---\n{res.stdout}\n"
            f"--- STDERR ---\n{res.stderr}",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(res.returncode)


def main():
    parser = argparse.ArgumentParser(
        description="Build GAFF2 topology from reference SDF and external atomic charges"
    )
    parser.add_argument(
        "--sdf", required=True, help="Input canonical reference 3D structure"
    )
    parser.add_argument(
        "--charges", required=True, help="Vector file of atomic charges (.crg)"
    )
    parser.add_argument("--charge", required=True, type=int, help="Net formal charge")
    parser.add_argument("--res_name", required=True, help="Residue name (3-4 chars)")
    parser.add_argument(
        "--forcefield", default="gaff2", help="Target force field (gaff or gaff2)"
    )
    parser.add_argument("--out_mol2", required=True, help="Output mol2 file")
    parser.add_argument("--out_frcmod", required=True, help="Output frcmod file")
    args = parser.parse_args()

    sdf_path = Path(args.sdf).resolve()
    charges_path = Path(args.charges).resolve()

    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Reference SDF not found or empty: {sdf_path}")
    if not charges_path.exists() or charges_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: External charges file not found or empty: {charges_path}"
        )

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(
            f"FATAL: Could not read valid molecular structure from {sdf_path}"
        )

    res_clean = args.res_name.strip().upper()[:4]
    if len(res_clean) < 3:
        res_clean = f"{res_clean:<3s}"

    # Ensure MDL V2000 format for complete Antechamber compatibility
    clean_v2000_sdf = "ligand_v2000.sdf"
    with Chem.SDWriter(clean_v2000_sdf) as writer:
        writer.SetForceV3000(False)
        writer.write(mol)

    run_cmd(
        [
            "antechamber",
            "-i",
            clean_v2000_sdf,
            "-fi",
            "sdf",
            "-o",
            args.out_mol2,
            "-fo",
            "mol2",
            "-at",
            args.forcefield,
            "-rn",
            res_clean,
            "-c",
            "rc",
            "-cf",
            str(charges_path),
            "-nc",
            str(args.charge),
            "-dr",
            "no",
        ],
        desc="Antechamber GAFF2 Atom Typing with External Charges",
    )

    run_cmd(
        [
            "parmchk2",
            "-i",
            args.out_mol2,
            "-f",
            "mol2",
            "-o",
            args.out_frcmod,
            "-s",
            "2",
        ],
        desc="Parmchk2 Parameter File Generation",
    )

    out_mol2_path = Path(args.out_mol2)
    out_frcmod_path = Path(args.out_frcmod)

    if not out_mol2_path.exists() or out_mol2_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Antechamber failed to produce output mol2: {args.out_mol2}"
        )
    if not out_frcmod_path.exists() or out_frcmod_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Parmchk2 failed to produce output frcmod: {args.out_frcmod}"
        )

    print(
        f"✓ Successfully generated GAFF2 topology: {out_mol2_path.name} and {out_frcmod_path.name}",
        flush=True,
    )


if __name__ == "__main__":
    main()
