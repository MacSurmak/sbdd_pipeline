#!/usr/bin/env python3
"""
bin/merge_complex.py
Cleanly merges receptor PDB with organic ligands and cofactors.
Ensures distinct residue numbers, unique chain IDs, and proper TER cards for AMBER tLEaP.
"""
import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Merge receptor with organics for Amber"
    )
    parser.add_argument("--rec_pdb", required=True, help="Input receptor PDB")
    parser.add_argument(
        "--lig_sdf", help="Active ligand SDF file (if holo)"
    )
    parser.add_argument("--out_pdb", required=True, help="Output complex PDB")
    args = parser.parse_args()

    rec_path = Path(args.rec_pdb)
    if not rec_path.exists():
        print(f"FATAL: Receptor PDB {rec_path} not found.", file=sys.stderr)
        sys.exit(1)

    rec_lines = rec_path.read_text(errors="ignore").splitlines()
    max_res_seq = 0
    used_chains = set()
    out_lines = []

    # 1. Process receptor atoms & keep original TER cards
    for line in rec_lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                res_seq = int(line[22:26].strip())
                if res_seq > max_res_seq:
                    max_res_seq = res_seq
                chain = line[21]
                if chain.strip():
                    used_chains.add(chain)
            except ValueError:
                pass
            out_lines.append(line)
        elif line.startswith("TER"):
            out_lines.append(line)

    # Ensure a TER separates receptor from any downstream organics
    if out_lines and not out_lines[-1].startswith("TER"):
        out_lines.append("TER")

    # Available chain IDs for small molecules (working backwards from Z)
    available_chains = [c for c in "ZYXWVUTSRQPONMLKJIHGFEDCB" if c not in used_chains]

    # 2. Append active ligand as an independent, TER-separated residue
    if args.lig_sdf:
        p_sdf = Path(args.lig_sdf)
        if not p_sdf.exists() or p_sdf.stat().st_size == 0:
            print(f"FATAL: Ligand SDF is missing or empty: {p_sdf}", file=sys.stderr)
            sys.exit(1)
        
        try:
            from rdkit import Chem
            suppl = Chem.SDMolSupplier(str(p_sdf), removeHs=False)
            mol = next((m for m in suppl if m is not None), None)
            if not mol:
                raise ValueError("RDKit failed to parse ligand.")
            
            lig_pdb_block = Chem.MolToPDBBlock(mol)
            
            max_res_seq += 1
            chain_id = available_chains.pop(0) if available_chains else "Z"

            for line in lig_pdb_block.splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    # Ensure correct PDB columns:
                    # 0-20: Record name, serial, atom name, altLoc, resName
                    # 21: Chain ID
                    # 22-25: resSeq
                    # 26+: iCode, coords, occupancy, tempFactor, element
                    atom_part = line[:21]
                    rest = line[26:] if len(line) >= 26 else ""
                    new_line = f"{atom_part[:21]:<21}{chain_id}{max_res_seq:4d}{rest}"
                    out_lines.append(new_line)
            
            out_lines.append("TER")
        except Exception as e:
            print(f"FATAL: Could not convert ligand SDF to PDB block: {e}", file=sys.stderr)
            sys.exit(1)

    out_lines.append("END")
    Path(args.out_pdb).write_text("\n".join(out_lines) + "\n")
    lig_msg = "with 1 active ligand." if args.lig_sdf else "(Apo)."
    print(f"Successfully assembled {args.out_pdb} {lig_msg}")


if __name__ == "__main__":
    main()
