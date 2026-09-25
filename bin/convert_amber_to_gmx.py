#!/usr/bin/env python3
"""
bin/convert_amber_to_gmx.py
Converts Amber PRMTOP and INPCRD files to GROMACS topology and GRO file using ParmEd.
"""
import argparse
import sys


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

    print("Saving GROMACS step3_input.gro...")
    struct.save("step3_input.gro", overwrite=True)
    print("Conversion successful.")


if __name__ == "__main__":
    main()
