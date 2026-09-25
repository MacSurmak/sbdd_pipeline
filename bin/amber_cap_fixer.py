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
        import numpy as np
    except ImportError:
        print(
            "ERROR: ParmEd and numpy are required but not installed.", file=sys.stderr
        )
        sys.exit(1)

    struct = pmd.load_file(args.input_pdb)

    # 1. Identify Disulfide Bridges (SG-SG distance < 2.5 A)
    cyx_res = set()
    cys_residues = [r for r in struct.residues if r.name in ["CYS", "CYX"]]
    for i, res1 in enumerate(cys_residues):
        sg1 = next((a for a in res1.atoms if a.name == "SG"), None)
        if not sg1:
            continue
        for res2 in cys_residues[i + 1 :]:
            sg2 = next((a for a in res2.atoms if a.name == "SG"), None)
            if not sg2:
                continue
            dist = np.linalg.norm(
                np.array([sg1.xx, sg1.xy, sg1.xz]) - np.array([sg2.xx, sg2.xy, sg2.xz])
            )
            if dist < 2.5:
                cyx_res.add(res1)
                cyx_res.add(res2)

    for res in struct.residues:
        if res.name == "ACE":
            for a in res.atoms:
                if a.name in ["1H", "1HA", "H1", "HH31"]:
                    a.name = "H1"
                elif a.name in ["2H", "2HA", "H2", "HH32"]:
                    a.name = "H2"
                elif a.name in ["3H", "3HA", "H3", "HH33"]:
                    a.name = "H3"
                elif a.name == "CA":
                    a.name = "CH3"
        elif res.name in ["NMA", "NME"]:
            res.name = "NME"
            for a in res.atoms:
                if a.name in ["CA", "CH3"]:
                    a.name = "C"
                elif a.name in ["1HA", "1H", "H1", "HH31"]:
                    a.name = "H1"
                elif a.name in ["2HA", "2H", "H2", "HH32"]:
                    a.name = "H2"
                elif a.name in ["3HA", "3H", "H3", "HH33"]:
                    a.name = "H3"
                elif a.name == "HN":
                    a.name = "H"
        elif res.name == "ASP":
            if any(a.name in ["HD2", "2HD"] for a in res.atoms):
                res.name = "ASH"
        elif res.name == "GLU":
            if any(a.name in ["HE2", "2HE"] for a in res.atoms):
                res.name = "GLH"
        elif res.name == "LYS":
            if not any(a.name in ["HZ3", "3HZ", "HZ"] for a in res.atoms):
                res.name = "LYN"
        elif res.name in ["HIS", "HSD", "HSE", "HSP"]:
            names = [a.name for a in res.atoms]
            if "HD1" in names and "HE2" in names:
                res.name = "HIP"
            elif "HD1" in names:
                res.name = "HID"
            elif "HE2" in names:
                res.name = "HIE"
            else:
                res.name = "HIE"
        elif res in cyx_res:
            res.name = "CYX"
            hg_atom = next((a for a in res.atoms if a.name == "HG"), None)
            if hg_atom:
                struct.atoms.remove(hg_atom)
                res.atoms.remove(hg_atom)
        else:
            # Fix standard backbone amide protons from Schrodinger format
            has_2h = any(a.name in ["2H", "H2"] for a in res.atoms)
            has_3h = any(a.name in ["3H", "H3"] for a in res.atoms)
            for a in res.atoms:
                if a.name == "HN":
                    a.name = "H"
                elif a.name == "1H":
                    if has_2h or has_3h:
                        a.name = "H1"
                    else:
                        a.name = "H"
                elif a.name == "2H":
                    a.name = "H2"
                elif a.name == "3H":
                    a.name = "H3"

    struct.save(args.output_pdb, overwrite=True)
    print(f"Successfully standardized terminal caps and saved to {args.output_pdb}")


if __name__ == "__main__":
    main()
