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
        "--lig_sdf",
        default=None,
        help="Active ligand SDF file (strictly for holo systems)",
    )
    parser.add_argument(
        "--lig_resname",
        default=None,
        help="Strict PDB residue name for ligand (strictly for holo systems)",
    )
    parser.add_argument(
        "--ion_sdfs",
        nargs="*",
        default=[],
        help="SDF files strictly representing ions (from 3NF DB)",
    )
    parser.add_argument(
        "--ion_resnames",
        nargs="*",
        default=[],
        help="Residue names corresponding to ion SDFs",
    )
    parser.add_argument(
        "--cofactor_sdfs",
        nargs="*",
        default=[],
        help="SDF files representing organic cofactors",
    )
    parser.add_argument(
        "--cofactor_resnames",
        nargs="*",
        default=[],
        help="Residue names corresponding to cofactor SDFs",
    )
    parser.add_argument("--out_pdb", required=True, help="Output complex PDB")
    args = parser.parse_args()

    # Strict Fail-Fast: mutually conditional validation between ligand coordinates and residue name
    if args.lig_sdf and not args.lig_resname:
        print(
            "FATAL: --lig_resname is strictly required when --lig_sdf is provided for a holo system!",
            file=sys.stderr,
        )
        sys.exit(1)
    if not args.lig_sdf and args.lig_resname:
        print(
            "FATAL: --lig_resname cannot be specified for an apo system without --lig_sdf!",
            file=sys.stderr,
        )
        sys.exit(1)

    if len(args.ion_sdfs) != len(args.ion_resnames):
        print(
            f"FATAL: Mismatch between --ion_sdfs ({len(args.ion_sdfs)}) and --ion_resnames ({len(args.ion_resnames)})!",
            file=sys.stderr,
        )
        sys.exit(1)

    if len(args.cofactor_sdfs) != len(args.cofactor_resnames):
        print(
            f"FATAL: Mismatch between --cofactor_sdfs ({len(args.cofactor_sdfs)}) and --cofactor_resnames ({len(args.cofactor_resnames)})!",
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
    topology_lines = []

    skip_resnames = set()
    if args.lig_resname:
        skip_resnames.add(args.lig_resname.strip().upper())
    for r in args.ion_resnames:
        skip_resnames.add(r.strip().upper())
    for r in args.cofactor_resnames:
        skip_resnames.add(r.strip().upper())

    # 1. Process receptor atoms, filter out existing cofactors/ligands, and preserve covalent topology mapping
    for line in rec_lines:
        if line.startswith(("ATOM", "HETATM")):
            r_name = line[17:20].strip().upper()
            if r_name in skip_resnames:
                continue

            try:
                res_seq = int(line[22:26].strip())
            except ValueError as e:
                print(
                    f"FATAL: Malformed residue sequence number in line '{line}': {e}",
                    file=sys.stderr,
                )
                sys.exit(1)
            if res_seq > max_res_seq:
                max_res_seq = res_seq
            chain = line[21]
            if chain.strip():
                used_chains.add(chain)
            out_lines.append(line)
        elif line.startswith("TER"):
            out_lines.append(line)
        elif line.startswith(("CONECT", "LINK", "SSBOND")):
            topology_lines.append(line)

    # Ensure a TER separates receptor from any downstream organics
    if out_lines and not out_lines[-1].startswith("TER"):
        out_lines.append("TER")

    # Available chain IDs for small molecules (working backwards from Z)
    available_chains = [c for c in "ZYXWVUTSRQPONMLKJIHGFEDCB" if c not in used_chains]

    from rdkit import Chem

    # 2. Helper to append small molecule residues with unique chains and TER cards
    def _append_small_molecule(sdf_path_str: str, resname: str, is_ion: bool = False):
        nonlocal max_res_seq
        p_sdf = Path(sdf_path_str)
        if not p_sdf.exists() or p_sdf.stat().st_size == 0:
            print(
                f"FATAL: Small molecule SDF missing or empty: {p_sdf}", file=sys.stderr
            )
            sys.exit(1)

        suppl = Chem.SDMolSupplier(str(p_sdf), removeHs=False, sanitize=False)
        mol = next((m for m in suppl if m is not None), None)
        if not mol:
            print(
                f"FATAL: RDKit failed to parse molecule from {p_sdf}", file=sys.stderr
            )
            sys.exit(1)
        mol.UpdatePropertyCache(strict=False)

        if not available_chains:
            print(
                "FATAL: Exhausted all available chain IDs for small molecules!",
                file=sys.stderr,
            )
            sys.exit(1)

        chain_id = available_chains.pop(0)
        max_res_seq += 1

        res_code = resname.strip().upper()[:4]
        if len(res_code) < 3:
            res_code = f"{res_code:<3s}"

        if is_ion:
            # Deterministic, non-heuristic branch: entity is strictly an ion from SSOT DB
            if mol.GetNumAtoms() != 1:
                print(
                    f"FATAL: Ion entity '{resname}' in {p_sdf.name} contains {mol.GetNumAtoms()} atoms! "
                    f"Monoatomic ions in AMBER must contain strictly 1 atom.",
                    file=sys.stderr,
                )
                sys.exit(1)

            conf = mol.GetConformer()
            pos = conf.GetAtomPosition(0)
            elem = mol.GetAtomWithIdx(0).GetSymbol().upper()
            clean_atom = res_code.strip()
            # Standard PDB column 13-16 alignment strictly matching user-declared ion name
            formatted_atom = (
                f" {clean_atom:<3s}" if len(clean_atom) < 4 else f"{clean_atom:<4s}"
            )
            new_line = (
                f"HETATM{len(out_lines)+1:5d} {formatted_atom:4s} {res_code:>3s} "
                f"{chain_id}{max_res_seq:4d}    {pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00          {elem:>2s}"
            )
            out_lines.append(new_line)
        else:
            # Organic ligand or cofactor: preserve RDKit atom names matching Antechamber GAFF2 mol2
            pdb_block = Chem.MolToPDBBlock(mol)
            for line in pdb_block.splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    atom_name = line[12:16]
                    coords_and_rest = line[26:] if len(line) >= 26 else ""
                    new_line = f"HETATM{len(out_lines)+1:5d} {atom_name:<4s} {res_code:>3s} {chain_id}{max_res_seq:4d}{coords_and_rest}"
                    out_lines.append(new_line)

        out_lines.append("TER")

    # Append active ligand if present (strictly organic)
    if args.lig_sdf:
        _append_small_molecule(args.lig_sdf, args.lig_resname, is_ion=False)

    # Append ions (strictly typed from DB SSOT)
    for i_sdf, i_res in zip(args.ion_sdfs, args.ion_resnames):
        _append_small_molecule(i_sdf, i_res, is_ion=True)

    # Append organic cofactors (strictly typed from DB SSOT)
    for c_sdf, c_res in zip(args.cofactor_sdfs, args.cofactor_resnames):
        _append_small_molecule(c_sdf, c_res, is_ion=False)

    out_lines.extend(topology_lines)
    out_lines.append("END")
    Path(args.out_pdb).write_text("\n".join(out_lines) + "\n")
    print(
        f"Successfully assembled {args.out_pdb} with {bool(args.lig_sdf)} ligand, {len(args.ion_sdfs)} ion(s), and {len(args.cofactor_sdfs)} cofactor(s)."
    )


if __name__ == "__main__":
    main()
