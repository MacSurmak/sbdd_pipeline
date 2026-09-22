#!/usr/bin/env python3
"""
File: src/utils/merge_complex_parmed.py
Description: Safely merges Protein PDB and Ligand MOL2 using ParmEd.
Converts ResidueTemplate to Structure if necessary.
"""
import sys
import parmed as pmd
from loguru import logger

def merge_complex(prot_pdb: str, lig_mol2: str, out_pdb: str):
    logger.info(f"Merging {prot_pdb} and {lig_mol2}...")
    prot = pmd.load_file(prot_pdb)
    lig = pmd.load_file(lig_mol2)

    # Convert ResidueTemplate to Structure to expose .residues
    if hasattr(lig, "to_structure"):
        lig = lig.to_structure()

    # Find a free chain ID for the ligand
    used_chains = {res.chain for res in prot.residues}
    lig_chain = next((c for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if c not in used_chains), "Z")

    for res in lig.residues:
        res.chain = lig_chain
        res.name = "UNK"

    complex_struct = prot + lig
    complex_struct.save(out_pdb, overwrite=True)
    logger.info(f"Complex saved to {out_pdb} (Ligand assigned to chain {lig_chain})")

if __name__ == "__main__":
    merge_complex(sys.argv[1], sys.argv[2], sys.argv[3])