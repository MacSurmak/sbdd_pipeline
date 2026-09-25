#!/usr/bin/env python3
"""
bin/amber_cap_fixer.py
Uses official ParmEd API to translate Schrödinger terminal cap nomenclature (ACE/NMA) 
into AMBER standard nomenclature (ACE/NME) so tLEaP can build topology.
Strictly preserves all coordinates, geometries, and protonation states.
"""
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Standardize terminal caps for AMBER")
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("output_pdb", help="Output PDB file")
    args = parser.parse_args()

    try:
        import parmed as pmd
    except ImportError:
        print("ERROR: ParmEd is required but not installed.", file=sys.stderr)
        sys.exit(1)

    struct = pmd.load_file(args.input_pdb)

    for res in struct.residues:
        if res.name == "ACE":
            for a in res.atoms:
                if a.name in ["1H", "1HA", "H1"]:
                    a.name = "H1"
                elif a.name in ["2H", "2HA", "H2"]:
                    a.name = "H2"
                elif a.name in ["3H", "3HA", "H3"]:
                    a.name = "H3"
                elif a.name == "CA":
                    a.name = "CH3"
        elif res.name in ["NMA", "NME"]:
            res.name = "NME"
            for a in res.atoms:
                if a.name in ["CA", "CH3"]:
                    a.name = "C"
                elif a.name in ["1HA", "1H"]:
                    a.name = "H1"
                elif a.name in ["2HA", "2H"]:
                    a.name = "H2"
                elif a.name in ["3HA", "3H"]:
                    a.name = "H3"
                elif a.name == "HN":
                    a.name = "H"
        else:
            # Fix standard backbone amide protons from Schrodinger format
            for a in res.atoms:
                if a.name in ["HN", "1H"]:
                    a.name = "H"

    struct.save(args.output_pdb, overwrite=True)
    print(f"Successfully standardized terminal caps and saved to {args.output_pdb}")


if __name__ == "__main__":
    main()
