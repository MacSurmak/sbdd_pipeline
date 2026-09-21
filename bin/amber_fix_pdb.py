#!/usr/bin/env python3
"""
File: src/utils/amber_prep_parmed.py
Description: Robust preparation of Schrodinger PDBs for AMBER using ParmEd.
Handles NMA->NME, ACE hydrogens, HIS protonation, and CYS->CYX mapping.
"""
import sys
import numpy as np
import parmed as pmd
from loguru import logger

def fix_for_amber(in_pdb: str, out_pdb: str):
    logger.info(f"Loading {in_pdb} via ParmEd...")
    struct = pmd.load_file(in_pdb)

    # 1. Identify Disulfide Bridges (SG-SG distance < 2.5 A)
    cyx_res = set()
    cys_residues = [r for r in struct.residues if r.name in ["CYS", "CYX"]]
    for i, res1 in enumerate(cys_residues):
        sg1 = next((a for a in res1.atoms if a.name == "SG"), None)
        if not sg1: continue
        for res2 in cys_residues[i+1:]:
            sg2 = next((a for a in res2.atoms if a.name == "SG"), None)
            if not sg2: continue
            dist = np.linalg.norm(np.array([sg1.xx, sg1.xy, sg1.xz]) - np.array([sg2.xx, sg2.xy, sg2.xz]))
            if dist < 2.5:
                cyx_res.add(res1)
                cyx_res.add(res2)

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
                if a.name in ["CA", "CH3"]: a.name = "C"
                elif a.name in ["1HA", "1H", "H1", "HH31"]: a.name = "H1"
                elif a.name in ["2HA", "2H", "H2", "HH32"]: a.name = "H2"
                elif a.name in ["3HA", "3H", "H3", "HH33"]: a.name = "H3"
        
        # N-Terminal Cap (ACE: atoms C, O, CH3, H1, H2, H3)
        elif res.name == "ACE":
            for a in res.atoms:
                if a.name in ["1H", "H1", "HH31"]: a.name = "H1"
                elif a.name in ["2H", "H2", "HH32"]: a.name = "H2"
                elif a.name in ["3H", "H3", "HH33"]: a.name = "H3"
        
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
            if "HD1" in names and "HE2" in names: res.name = "HIP"
            elif "HD1" in names: res.name = "HID"
            elif "HE2" in names: res.name = "HIE"
            else: res.name = "HIE"
        
        # Cysteine Disulfides
        elif res in cyx_res:
            res.name = "CYX"
            hg_atom = next((a for a in res.atoms if a.name == "HG"), None)
            if hg_atom:
                struct.atoms.remove(hg_atom)
                res.atoms.remove(hg_atom)

    logger.info(f"Saving AMBER-compliant PDB to {out_pdb}...")
    struct.save(out_pdb, overwrite=True)

if __name__ == "__main__":
    fix_for_amber(sys.argv[1], sys.argv[2])