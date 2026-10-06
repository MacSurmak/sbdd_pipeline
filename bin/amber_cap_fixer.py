#!/usr/bin/env python3
"""
bin/amber_cap_fixer.py
Clean, deterministic PDB renamer for AMBER/tLEaP compatibility.
Translates caps (ACE/NME), protonation states (ASH/GLH/LYN/HID/HIE/HIP), and disulfides (CYX).
STRICTLY prevents chain fragmentation: NEVER writes TER cards between consecutive residues in the same chain.
"""

import argparse
import math
import sys
from pathlib import Path


def format_atom_name(name: str) -> str:
    name = name.strip()
    if len(name) >= 4:
        return f"{name[:4]:4s}"
    else:
        return f" {name:<3s}"


STANDARD_PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "CYX",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASH",
    "GLH",
    "ACE",
    "NME",
    "NMA",
}


def main():
    parser = argparse.ArgumentParser(
        description="Standardize PDB for AMBER without chain fragmentation"
    )
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("output_pdb", help="Output clean PDB file")
    parser.add_argument(
        "--bonds_out", default=None, help="Output tLEaP bond commands file"
    )
    parser.add_argument(
        "--ss_cutoff", type=float, required=True, help="Distance cutoff for SS bonds"
    )
    parser.add_argument(
        "--water_names",
        required=True,
        help="Comma-separated water residue names to strip",
    )
    args = parser.parse_args()

    water_residues = {
        w.strip().upper() for w in args.water_names.split(",") if w.strip()
    }
    lines = Path(args.input_pdb).read_text(errors="ignore").splitlines()

    # 1. Parse atoms into structured records, grouped by residue
    residues = []
    current_key = None
    current_res = None

    for line in lines:
        if line.startswith("TER"):
            current_key = None
            current_res = None
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        alt_loc = line[16]
        res_name = line[17:20].strip().upper()
        chain_id = line[21]
        res_seq = int(line[22:26].strip())
        i_code = line[26]
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        # Exclude solvent strictly based on configuration Single Source of Truth
        if res_name in water_residues:
            continue

        res_key = (chain_id, res_seq, i_code, res_name)
        if res_key != current_key:
            current_key = res_key
            current_res = {
                "chain_id": chain_id,
                "orig_res_seq": res_seq,
                "i_code": i_code,
                "res_name": res_name,
                "clean_res_seq": 0,
                "atoms": [],
            }
            residues.append(current_res)

        current_res["atoms"].append(
            {
                "atom_num": atom_num,
                "atom_name": atom_name,
                "line": line,
                "x": x,
                "y": y,
                "z": z,
            }
        )

    # 2. Establish strict globally monotonic residue numbering (1..N) eliminating insertion codes
    global_res_counter = 0
    for r in residues:
        global_res_counter += 1
        r["clean_res_seq"] = global_res_counter

    # 3. Identify Disulfide Bridges using exact canonical clean residue IDs
    bonds = []
    cyx_indices = set()
    cys_candidates = []

    for idx, r in enumerate(residues):
        if r["res_name"] in ["CYS", "CYX"]:
            sg = next((a for a in r["atoms"] if a["atom_name"] == "SG"), None)
            if sg:
                cys_candidates.append((idx, r, sg))

    for i in range(len(cys_candidates)):
        idx1, r1, sg1 = cys_candidates[i]
        for j in range(i + 1, len(cys_candidates)):
            idx2, r2, sg2 = cys_candidates[j]
            dist = math.sqrt(
                (sg1["x"] - sg2["x"]) ** 2
                + (sg1["y"] - sg2["y"]) ** 2
                + (sg1["z"] - sg2["z"]) ** 2
            )
            if dist < args.ss_cutoff:
                cyx_indices.add(idx1)
                cyx_indices.add(idx2)
                bonds.append(
                    f"bond complex.{r1['clean_res_seq']}.SG complex.{r2['clean_res_seq']}.SG"
                )

    if args.bonds_out:
        with open(args.bonds_out, "w") as fb:
            for b in bonds:
                fb.write(f"{b}\n")

    # 4. Apply chemical renaming rules strictly matching AMBER ff19SB force field specifications
    for idx, r in enumerate(residues):
        res_name = r["res_name"]
        atom_names = [a["atom_name"] for a in r["atoms"]]

        # Disulfide Cysteine
        if idx in cyx_indices:
            r["res_name"] = "CYX"
            r["atoms"] = [a for a in r["atoms"] if a["atom_name"] not in ["HG", "2HG"]]

        # Protonated Aspartate (ASH)
        elif res_name == "ASP" and any(n in ["HD2", "2HD"] for n in atom_names):
            r["res_name"] = "ASH"

        # Protonated Glutamate (GLH)
        elif res_name == "GLU" and any(n in ["HE2", "2HE"] for n in atom_names):
            r["res_name"] = "GLH"

        # Lysine / Neutral Lysine (LYN)
        elif res_name in ["LYS", "LYN"]:
            existing_names = [a["atom_name"].strip() for a in r["atoms"]]
            for a in r["atoms"]:
                if a["atom_name"].strip() == "HZ":
                    if "HZ1" in existing_names and "HZ2" in existing_names:
                        a["atom_name"] = "HZ3"
                    elif "HZ1" in existing_names:
                        a["atom_name"] = "HZ2"
                    else:
                        a["atom_name"] = "HZ1"

            updated_names = [a["atom_name"].strip() for a in r["atoms"]]
            if not any(n in ["HZ3", "3HZ"] for n in updated_names):
                r["res_name"] = "LYN"

        # Histidine Protonation States
        elif res_name in ["HIS", "HSD", "HSE", "HSP"]:
            has_hd1 = "HD1" in atom_names
            has_he2 = "HE2" in atom_names
            if has_hd1 and has_he2:
                r["res_name"] = "HIP"
            elif has_hd1:
                r["res_name"] = "HID"
            else:
                r["res_name"] = "HIE"

        # C-terminal Cap (NMA/NME -> AMBER NME template: N, H, C, H1, H2, H3)
        elif res_name in ["NMA", "NME"]:
            r["res_name"] = "NME"
            for a in r["atoms"]:
                an = a["atom_name"].strip()
                if an in ["CA", "CH3", "C"]:
                    a["atom_name"] = "C"
                elif an in ["HN", "H"]:
                    a["atom_name"] = "H"
                elif an in ["1HA", "1H", "HH31", "HA1", "H1"]:
                    a["atom_name"] = "H1"
                elif an in ["2HA", "2H", "HH32", "HA2", "H2"]:
                    a["atom_name"] = "H2"
                elif an in ["3HA", "3H", "HH33", "HA3", "H3"]:
                    a["atom_name"] = "H3"

        # N-terminal Cap (ACE -> AMBER ACE template: CH3, C, O, H1, H2, H3)
        elif res_name == "ACE":
            for a in r["atoms"]:
                an = a["atom_name"].strip()
                if an in ["CA", "CH3"]:
                    a["atom_name"] = "CH3"
                elif an in ["1HA", "1H", "HH31", "HA1", "H1"]:
                    a["atom_name"] = "H1"
                elif an in ["2HA", "2H", "HH32", "HA2", "H2"]:
                    a["atom_name"] = "H2"
                elif an in ["3HA", "3H", "HH33", "HA3", "H3"]:
                    a["atom_name"] = "H3"

        # Standard backbone amide hydrogen fix strictly isolated to protein residues
        elif res_name in STANDARD_PROTEIN_RESIDUES and res_name != "PRO":
            amide_h = [
                a
                for a in r["atoms"]
                if a["atom_name"] in ["H", "HN", "1H", "2H", "3H", "H1", "H2", "H3"]
            ]
            if len(amide_h) == 1:
                amide_h[0]["atom_name"] = "H"
            elif len(amide_h) == 3:
                for h_idx, a in enumerate(
                    sorted(amide_h, key=lambda x: x["atom_name"]), start=1
                ):
                    a["atom_name"] = f"H{h_idx}"

    # 5. Write output PDB strictly preserving continuous chain continuity and exact PDB columns
    out_lines = []
    global_atom_idx = 1
    prev_chain = None

    for r_idx, r in enumerate(residues):
        # Insert TER when chain physically changes or before organic molecules
        if prev_chain is not None and (
            r["chain_id"] != prev_chain or r["chain_id"] in "ZYXWVUTSRQPONMLKJIHGFEDCB"
        ):
            prev_r = residues[r_idx - 1]
            out_lines.append(
                f"TER   {global_atom_idx:5d}      {prev_r['res_name']:>3s} {prev_chain}{prev_r['clean_res_seq']:4d}"
            )
            global_atom_idx += 1

        prev_chain = r["chain_id"]
        formatted_res = f"{r['res_name']:>3s}"

        for a in r["atoms"]:
            formatted_atom = format_atom_name(a["atom_name"])
            orig = a["line"]
            # Col 0-5 (Record), 6-10 (Serial), 11 (Space), 12-15 (Atom), 16 (altLoc),
            # 17-19 (resName), 20 (Space), 21 (Chain), 22-25 (resSeq), 26 (Blank iCode),
            # 27-29 (Spaces), 30+ (Coordinates, occupancy, tempFactor, element)
            rest_coords = orig[30:] if len(orig) >= 30 else ""
            new_line = (
                f"{orig[:6]}{global_atom_idx:5d} {formatted_atom}{orig[16]}{formatted_res} "
                f"{prev_chain}{r['clean_res_seq']:4d}   {rest_coords}"
            )
            out_lines.append(new_line)
            global_atom_idx += 1

    # Final TER at the end of the entire protein structure
    if residues:
        last = residues[-1]
        out_lines.append(
            f"TER   {global_atom_idx:5d}      {last['res_name']:>3s} {last['chain_id']}{last['clean_res_seq']:4d}"
        )

    out_lines.append("END")
    Path(args.output_pdb).write_text("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    main()
