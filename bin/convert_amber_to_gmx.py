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
    parser.add_argument("--inpcrd", required=True, help="Input Amber inpcrd file")
    args = parser.parse_args()

    try:
        import parmed as pmd
    except ImportError:
        print(
            "Error: ParmEd is not installed in the current environment.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Loading {args.prmtop} and {args.inpcrd}...")
    struct = pmd.load_file(args.prmtop, args.inpcrd)

    print("Saving GROMACS topol.top...")
    struct.save("topol.top", format="gromacs", overwrite=True)

    # Fail-Fast Validation: Ensure terminal caps are covalently integrated and did not fragment
    topol_text = Path("topol.top").read_text(encoding="utf-8", errors="ignore")
    isolated_caps = []
    for match in re.finditer(
        r"\[\s*moleculetype\s*\]\s*\n\s*;[^\n]*\n\s*(\S+)", topol_text
    ):
        mname = match.group(1).strip()
        if mname in ["NME", "NMA", "ACE"]:
            isolated_caps.append(mname)

    if isolated_caps:
        print(
            f"FATAL: TopologyIntegrityError: Terminal cap(s) {isolated_caps} were exported as detached, "
            f"independent moleculetype(s) in topol.top! Peptide bond to protein chain was not formed in tLEaP. "
            f"Execution halted by Fail-Fast law.",
            file=sys.stderr,
        )
        sys.exit(1)

    print("Saving GROMACS step3_input.gro...")
    struct.save("step3_input.gro", overwrite=True)
    print("Conversion successful.")


if __name__ == "__main__":
    main()
