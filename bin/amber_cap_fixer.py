#!/usr/bin/env python3
"""
bin/amber_cap_fixer.py
Uses official ParmEd API to translate Schrödinger terminal cap nomenclature (ACE/NMA) 
into AMBER standard nomenclature (ACE/NME) so tLEaP can build topology.
Strictly preserves all coordinates, geometries, and protonation states.
Restored to the robust logic of the original SBDD pipeline.
"""
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Standardize terminal caps for AMBER")
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("output_pdb", help="Output PDB file")
    parser.add_argument(
        "--bonds_out", default=None, help="Output tLEaP bond commands file"
    )
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
    bonds = []
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
                bonds.append(f"bond complex.{res1.number}.SG complex.{res2.number}.SG")

    if args.bonds_out:
        with open(args.bonds_out, "w") as f_bonds:
            for b in bonds:
                f_bonds.write(f"{b}\n")

    # 2. Rename residues and atoms
    for res in struct.residues:
        # Standard amino acid backbone amide hydrogen fix (H1 / HN / 1H -> H)
        if res.name not in ["ACE", "NME", "NMA"]:
            for a in res.atoms:
                if a.name in ["H1", "HN", "1H"]:
                    a.name = "H"

        # C-Terminal Cap (NMA / NME -> NME: atoms N, H, CH3, H1, H2, H3)
        if res.name in ["NMA", "NME"]:
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

        # N-Terminal Cap (ACE: atoms C, O, CH3, H1, H2, H3)
        elif res.name == "ACE":
            for a in res.atoms:
                if a.name in ["1H", "1HA", "H1", "HH31"]:
                    a.name = "H1"
                elif a.name in ["2H", "2HA", "H2", "HH32"]:
                    a.name = "H2"
                elif a.name in ["3H", "3HA", "H3", "HH33"]:
                    a.name = "H3"
                elif a.name == "CA":
                    a.name = "CH3"

        # Protonated Aspartate (ASP with HD2 -> ASH)
        elif res.name == "ASP":
            if any(a.name in ["HD2", "2HD"] for a in res.atoms):
                res.name = "ASH"

        # Protonated Glutamate (GLU with HE2 -> GLH)
        elif res.name == "GLU":
            if any(a.name in ["HE2", "2HE"] for a in res.atoms):
                res.name = "GLH"

        # Neutral Lysine (LYS without HZ3 -> LYN)
        elif res.name == "LYS":
            if not any(a.name in ["HZ3", "3HZ", "HZ"] for a in res.atoms):
                res.name = "LYN"

        # Histidine Protonation States (HIP, HID, HIE)
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

        # Cysteine Disulfides
        elif res in cyx_res:
            res.name = "CYX"
            hg_atom = next((a for a in res.atoms if a.name == "HG"), None)
            if hg_atom:
                struct.atoms.remove(hg_atom)
                res.atoms.remove(hg_atom)

    struct.save(args.output_pdb, overwrite=True)


if __name__ == "__main__":
    main()
