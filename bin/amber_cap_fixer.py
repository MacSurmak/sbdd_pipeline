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

    try:
        import parmed as pmd
    except ImportError:
        print("FATAL: ParmEd is required in current environment.", file=sys.stderr)
        sys.exit(1)

    # Sanitize legacy/obsolete PDB segids (columns 73-76) that cause ParmEd to split residues
    # when Maestro outputs heavy atoms with segid and added hydrogens with blank segid
    input_path = Path(args.input_pdb)
    raw_lines = input_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    has_segid = any(
        line.startswith(("ATOM", "HETATM")) and len(line) >= 76 and line[72:76].strip()
        for line in raw_lines
    )

    if has_segid:
        clean_lines = [
            (line[:72] + "    " + line[76:]) if (line.startswith(("ATOM", "HETATM")) and len(line) >= 76) else line
            for line in raw_lines
        ]
        temp_pdb = input_path.parent / f".tmp_{input_path.name}"
        temp_pdb.write_text("".join(clean_lines), encoding="utf-8")
        try:
            struct = pmd.load_file(str(temp_pdb))
        finally:
            temp_pdb.unlink(missing_ok=True)
    else:
        struct = pmd.load_file(args.input_pdb)

    # 1. Strip solvent strictly based on configuration SSOT via ParmEd masks
    water_residues = [
        w.strip().upper() for w in args.water_names.split(",") if w.strip()
    ]
    if water_residues:
        water_mask = ":" + ",".join(water_residues)
        struct.strip(water_mask)

    # 2. Strict globally monotonic residue numbering (1..N) eliminating insertion codes
    for idx, res in enumerate(struct.residues, start=1):
        res.number = idx
        res.insertion_code = ""

    # 3. Detect Disulfide Bridges geometrically via ParmEd atom positions
    cys_residues = [r for r in struct.residues if r.name in ("CYS", "CYX")]
    bonds = []
    cyx_residues = set()

    for i, r1 in enumerate(cys_residues):
        sg1 = next((a for a in r1.atoms if a.name == "SG"), None)
        if not sg1:
            continue
        for r2 in cys_residues[i + 1 :]:
            sg2 = next((a for a in r2.atoms if a.name == "SG"), None)
            if not sg2:
                continue
            dist = math.sqrt(
                (sg1.xx - sg2.xx) ** 2 + (sg1.xy - sg2.xy) ** 2 + (sg1.xz - sg2.xz) ** 2
            )
            if dist < args.ss_cutoff:
                cyx_residues.add(r1)
                cyx_residues.add(r2)
                bonds.append(f"bond complex.{r1.number}.SG complex.{r2.number}.SG")

    if args.bonds_out:
        with open(args.bonds_out, "w") as fb:
            for b in bonds:
                fb.write(f"{b}\n")

    # 4. Chemical Standardization: ParmEd object transformations matching ff19SB
    atoms_to_remove = []

    for idx, res in enumerate(struct.residues):
        # Ignore non-protein small molecules (ligands/cofactors): leave coordinates and names 100% untouched
        if res.name not in STANDARD_PROTEIN_RESIDUES and res not in cyx_residues:
            continue

        # Strictly enforce PDB v3.3 standard: all biopolymers and covalent caps MUST be ATOM records
        for a in res.atoms:
            a.hetatm = False

        # Ensure covalent caps are seamlessly integrated into the peptide chain (no TER cards)
        if res.name == "ACE":
            res.ter = False
            if idx + 1 < len(struct.residues):
                res.chain = struct.residues[idx + 1].chain

        if idx + 1 < len(struct.residues) and struct.residues[idx + 1].name in (
            "NMA",
            "NME",
        ):
            res.ter = False
            struct.residues[idx + 1].chain = res.chain
            # Strip terminal oxygen (OXT) if Maestro left it, as it chemically blocks C-N peptide bond formation
            oxt = next((a for a in res.atoms if a.name.strip() == "OXT"), None)
            if oxt:
                atoms_to_remove.append(oxt)

        atom_names = {a.name.strip() for a in res.atoms}

        # Disulfide Cysteine (remove thiol hydrogen)
        if res in cyx_residues:
            res.name = "CYX"
            for a in res.atoms:
                if a.name.strip() in ("HG", "2HG"):
                    atoms_to_remove.append(a)

        # Protonated Aspartate (ASH)
        elif res.name == "ASP" and ({"HD2", "2HD"} & atom_names):
            res.name = "ASH"

        # Protonated Glutamate (GLH)
        elif res.name == "GLU" and ({"HE2", "2HE"} & atom_names):
            res.name = "GLH"

        # Neutral Lysine (LYN)
        elif res.name in ("LYS", "LYN"):
            for a in res.atoms:
                if a.name.strip() == "HZ":
                    if "HZ1" in atom_names and "HZ2" in atom_names:
                        a.name = "HZ3"
                    elif "HZ1" in atom_names:
                        a.name = "HZ2"
                    else:
                        a.name = "HZ1"
            updated_names = {a.name.strip() for a in res.atoms}
            if not ({"HZ3", "3HZ"} & updated_names):
                res.name = "LYN"

        # Histidine Protonation States
        elif res.name in ("HIS", "HSD", "HSE", "HSP"):
            has_hd1 = "HD1" in atom_names
            has_he2 = "HE2" in atom_names
            if has_hd1 and has_he2:
                res.name = "HIP"
            elif has_hd1:
                res.name = "HID"
            else:
                res.name = "HIE"

        # C-terminal Cap (NMA/NME -> AMBER ff19SB template: N, H, C, H1, H2, H3)
        elif res.name in ("NMA", "NME"):
            res.name = "NME"
            for a in res.atoms:
                an = a.name.strip()
                if an in ("CA", "CH3", "C"):
                    a.name = "C"
                elif an in ("HN", "H"):
                    a.name = "H"
                elif an in ("1HA", "1H", "HH31", "HA1", "H1"):
                    a.name = "H1"
                elif an in ("2HA", "2H", "HH32", "HA2", "H2"):
                    a.name = "H2"
                elif an in ("3HA", "3H", "HH33", "HA3", "H3"):
                    a.name = "H3"

        # N-terminal Cap (ACE -> AMBER ff19SB template: CH3, C, O, H1, H2, H3)
        elif res.name == "ACE":
            for a in res.atoms:
                an = a.name.strip()
                if an in ("CA", "CH3"):
                    a.name = "CH3"
                elif an in ("1HA", "1H", "HH31", "HA1", "H1"):
                    a.name = "H1"
                elif an in ("2HA", "2H", "HH32", "HA2", "H2"):
                    a.name = "H2"
                elif an in ("3HA", "3H", "HH33", "HA3", "H3"):
                    a.name = "H3"

        # Backbone Amide Hydrogens (H, H1, H2, H3)
        elif res.name in STANDARD_PROTEIN_RESIDUES and res.name != "PRO":
            amide_h = [
                a
                for a in res.atoms
                if a.name.strip() in ("H", "HN", "1H", "2H", "3H", "H1", "H2", "H3")
            ]
            if len(amide_h) == 1:
                amide_h[0].name = "H"
            elif len(amide_h) == 3:
                for h_idx, a in enumerate(
                    sorted(amide_h, key=lambda x: x.name), start=1
                ):
                    a.name = f"H{h_idx}"

    # Prune detached hydrogens from CYX residues
    if atoms_to_remove:
        strip_indices = [str(a.idx + 1) for a in atoms_to_remove]
        struct.strip("@" + ",".join(strip_indices))

    # 5. Export clean PDB with strict chain and TER continuity via ParmEd C-engine
    struct.save(args.output_pdb, overwrite=True)


if __name__ == "__main__":
    main()
