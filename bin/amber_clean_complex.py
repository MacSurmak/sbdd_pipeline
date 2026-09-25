#!/usr/bin/env python3
"""
bin/amber_clean_complex.py
Cleans artifactual TER records from terminal caps (ACE/NME) and standardizes
protonated states (ASH, GLH) before feeding to Amber tLEaP.
Also detects disulfide bonds (SG-SG < 2.3A) and generates leap instructions.
"""
import argparse
import math
import sys


def main():
    parser = argparse.ArgumentParser(description="Clean complex PDB for Amber tLEaP")
    parser.add_argument("--complex_pdb", required=True, help="Input complex PDB")
    parser.add_argument(
        "--out_pdb", default="complex_clean.pdb", help="Output cleaned PDB"
    )
    parser.add_argument(
        "--out_bonds",
        default="disulfide_bonds.txt",
        help="Output tLEaP bond instructions",
    )
    args = parser.parse_args()

    try:
        with open(args.complex_pdb, "r") as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading {args.complex_pdb}: {e}", file=sys.stderr)
        sys.exit(1)

    cleaned = []
    last_chain = "A"
    cys_atoms = []

    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")) and len(l) >= 22:
            ch = l[21].strip()
            if ch:
                last_chain = ch
            if l[12:16].strip() == "SG" and l[17:20].strip() in ["CYS", "CYX", "CYM"]:
                cys_atoms.append(
                    {
                        "chain": l[21],
                        "resSeq": int(l[22:26]),
                        "x": float(l[30:38]),
                        "y": float(l[38:46]),
                        "z": float(l[46:54]),
                    }
                )

        # Remove artificial TER between ACE and NME caps
        if l.startswith("TER"):
            if i > 0 and "ACE" in lines[i - 1]:
                continue
            if i + 1 < len(lines) and "NME" in lines[i + 1]:
                continue

        # Inject correct chain ID to caps
        if ("ACE" in l or "NME" in l) and l.startswith(("ATOM", "HETATM")):
            l = l[:21] + last_chain + l[22:]

        # Standardize protonated Asp/Glu to Amber naming
        if "ASP" in l and ("HD2" in l or "2HD" in l):
            l = l[:17] + "ASH" + l[20:]
        elif "GLU" in l and ("HE2" in l or "2HE" in l):
            l = l[:17] + "GLH" + l[20:]

        cleaned.append(l)

    with open(args.out_pdb, "w") as f:
        f.writelines(cleaned)

    # Detect disulfide bonds within 2.3 Angstroms
    bonds = []
    for i in range(len(cys_atoms)):
        for j in range(i + 1, len(cys_atoms)):
            a, b = cys_atoms[i], cys_atoms[j]
            if a["resSeq"] == b["resSeq"] and a["chain"] == b["chain"]:
                continue
            dist = math.sqrt(
                (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2
            )
            if dist <= 2.3:
                bonds.append(f"bond complex.{a['resSeq']}.SG complex.{b['resSeq']}.SG")

    with open(args.out_bonds, "w") as fb:
        if bonds:
            fb.write("\n".join(bonds) + "\n")
        else:
            fb.write("\n")


if __name__ == "__main__":
    main()
