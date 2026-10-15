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
        "--lig_sdf", default=None, help="Active ligand SDF file (if holo)"
    )
    parser.add_argument(
        "--lig_resname", default="", help="Strict PDB residue name for ligand (if holo, e.g. KVR)"
    )
    parser.add_argument(
        "--extra_sdfs", nargs="*", default=[], help="Additional SDF files (cofactors, ions)"
    )
    parser.add_argument(
        "--extra_resnames", nargs="*", default=[], help="Residue names corresponding to extra SDFs"
    )
    parser.add_argument("--out_pdb", required=True, help="Output complex PDB")
    args = parser.parse_args()

    if len(args.extra_sdfs) != len(args.extra_resnames):
        print(
            f"FATAL: Mismatch between --extra_sdfs ({len(args.extra_sdfs)}) and --extra_resnames ({len(args.extra_resnames)})!",
            file=sys.stderr,
        )
        sys.exit(1)

    rec_path = Path(args.rec_pdb)
    if not rec_path.exists():
        print(f"FATAL: Receptor PDB {rec_path} not found.", file=sys.stderr)
        sys.exit(1)

    rec_lines = rec_path.read_text(errors="ignore").splitlines()
    max_res_seq = 0
    used_chains = set()
    out_lines = []

    skip_resnames = set()
    if args.lig_resname:
        skip_resnames.add(args.lig_resname.strip().upper())
    for r in args.extra_resnames:
        skip_resnames.add(r.strip().upper())

    # 1. Process receptor atoms, filter out cofactors/ligands already present, & keep original TER cards
    for line in rec_lines:
        if line.startswith(("ATOM", "HETATM")):
            r_name = line[17:20].strip().upper()
            if r_name in skip_resnames:
                continue

            try:
                res_seq = int(line[22:26].strip())
            except ValueError as e:
                print(f"FATAL: Malformed residue sequence number in line '{line}': {e}", file=sys.stderr)
                sys.exit(1)
            if res_seq > max_res_seq:
                max_res_seq = res_seq
            chain = line[21]
            if chain.strip():
                used_chains.add(chain)
            out_lines.append(line)
        elif line.startswith("TER"):
            out_lines.append(line)

    # Ensure a TER separates receptor from any downstream organics
    if out_lines and not out_lines[-1].startswith("TER"):
        out_lines.append("TER")

    # Available chain IDs for small molecules (working backwards from Z)
    available_chains = [c for c in "ZYXWVUTSRQPONMLKJIHGFEDCB" if c not in used_chains]

    from rdkit import Chem

    # 2. Helper to append small molecule residues with unique chains and TER cards
    def _append_small_molecule(sdf_path_str: str, resname: str):
        nonlocal max_res_seq
        p_sdf = Path(sdf_path_str)
        if not p_sdf.exists() or p_sdf.stat().st_size == 0:
            print(f"FATAL: Small molecule SDF missing or empty: {p_sdf}", file=sys.stderr)
            sys.exit(1)

        suppl = Chem.SDMolSupplier(str(p_sdf), removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if not mol:
            print(f"FATAL: RDKit failed to parse molecule from {p_sdf}", file=sys.stderr)
            sys.exit(1)

        pdb_block = Chem.MolToPDBBlock(mol)
        if not available_chains:
            print("FATAL: Exhausted all available chain IDs for small molecules!", file=sys.stderr)
            sys.exit(1)

        chain_id = available_chains.pop(0)
        max_res_seq += 1

        res_code = resname.strip().upper()[:4]
        if len(res_code) < 3:
            res_code = f"{res_code:<3s}"

        for line in pdb_block.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16]
                coords_and_rest = line[26:] if len(line) >= 26 else ""
                new_line = f"HETATM{len(out_lines)+1:5d} {atom_name:<4s} {res_code:>3s} {chain_id}{max_res_seq:4d}{coords_and_rest}"
                out_lines.append(new_line)

        out_lines.append("TER")

    # Append active ligand if present
    if args.lig_sdf:
        if not args.lig_resname:
            print("FATAL: --lig_resname must be provided when --lig_sdf is present.", file=sys.stderr)
            sys.exit(1)
        _append_small_molecule(args.lig_sdf, args.lig_resname)

    # Append all cofactors and ions
    for e_sdf, e_res in zip(args.extra_sdfs, args.extra_resnames):
        _append_small_molecule(e_sdf, e_res)

    out_lines.append("END")
    Path(args.out_pdb).write_text("\n".join(out_lines) + "\n")
    print(f"Successfully assembled {args.out_pdb} with {bool(args.lig_sdf)} ligand and {len(args.extra_sdfs)} cofactor(s)/ion(s).")


if __name__ == "__main__":
    main()
