#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
from colabfold.batch import convert_pdb_to_mmcif
from Bio.PDB import PDBParser, PDBIO

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", required=True)
    parser.add_argument("--fasta", required=True)
    args = parser.parse_args()

    os.makedirs("custom_templates", exist_ok=True)
    
    # Write fasta query
    with open(args.fasta, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith(">")]
    
    with open("query.fasta", "w") as f:
        f.write(">tmpl\n")
        f.write("".join(lines) + "\n")

    # Write template
    p = PDBParser(QUIET=True)
    st = p.get_structure('tmpl', args.protein)
    
    # Force chain ID to A for ColabFold mapping
    for model in st:
        for chain in model:
            chain.id = 'A'
            break
        break

    tmpl_pdb = Path("custom_templates/tmpl.pdb")
    io = PDBIO()
    io.set_structure(st)
    io.save(str(tmpl_pdb))

    convert_pdb_to_mmcif(tmpl_pdb)
    tmpl_pdb.unlink()

if __name__ == "__main__":
    main()
