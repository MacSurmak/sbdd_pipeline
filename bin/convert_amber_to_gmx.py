#!/usr/bin/env python3
"""
bin/convert_amber_to_gmx.py
Converts Amber PRMTOP and INPCRD files to GROMACS topology and GRO file using ParmEd.
"""
import argparse
import sys
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Convert Amber Topology to GROMACS via ParmEd"
    )
    parser.add_argument("--prmtop", required=True, help="Input Amber prmtop file")
    parser.add_argument("--inpcrd", default=None, help="Input Amber inpcrd file")
    parser.add_argument("--out_top", default="topol.top", help="Output GROMACS top file")
    parser.add_argument("--out_gro", default="step3_input.gro", help="Output GROMACS gro file")
    args = parser.parse_args()

    try:
        import parmed as pmd
    except ImportError:
        print(
            "Error: ParmEd is not installed in the current environment.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.inpcrd:
        print(f"Loading {args.prmtop} and {args.inpcrd}...")
        struct = pmd.load_file(args.prmtop, args.inpcrd)
    else:
        print(f"Loading dry topology {args.prmtop}...")
        struct = pmd.load_file(args.prmtop)

    print(f"Saving GROMACS {args.out_top}...")
    struct.save(args.out_top, format="gromacs", overwrite=True)

    # Fail-Fast Validation: Ensure terminal caps are covalently integrated and did not fragment
    topol_text = Path(args.out_top).read_text(encoding="utf-8", errors="ignore")
    isolated_caps = []
    in_moleculetype = False

    for raw_line in topol_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("["):
            in_moleculetype = "moleculetype" in line.lower()
            continue
        if in_moleculetype:
            if line.startswith(";") or line.startswith("#"):
                continue
            mname = line.split()[0].strip()
            if mname in ["NME", "NMA", "ACE"]:
                isolated_caps.append(mname)
            in_moleculetype = False

    if isolated_caps:
        print(
            f"FATAL: TopologyIntegrityError: Terminal cap(s) {isolated_caps} were exported as detached, "
            f"independent moleculetype(s) in topol.top! Peptide bond to protein chain was not formed in tLEaP. "
            f"Execution halted by Fail-Fast law.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.inpcrd:
        if struct.box is None or any(float(x) <= 0.0 for x in struct.box[:3]):
            print(
                f"FATAL: Periodic box vectors are missing or non-positive in {args.prmtop} / {args.inpcrd}: {struct.box}! "
                "Execution halted by Fail-Fast law.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(
            f"Saving GROMACS {args.out_gro} (box: {struct.box[0]:.2f} x {struct.box[1]:.2f} x {struct.box[2]:.2f} A)..."
        )
        struct.save(args.out_gro, overwrite=True)
    print("Conversion successful.")


if __name__ == "__main__":
    main()
