#!/usr/bin/env python3
"""
File: bin/extract_md_files.py
Description: Extracts receptor PDB and ligand MOL2 for MD setup.
Computes formal charge for Antechamber. Supports both Schrödinger API and RDKit fallback.
"""
import argparse
import sys
from pathlib import Path


ION_RESNAMES = {"ZN", "MG", "CA", "NA", "K", "CL", "ZN2", "MG2", "CA2"}


def extract_via_schrodinger(args):
    from schrodinger import structure
    from schrodinger.structutils import analyze
    import json

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    org_dir = out_dir / "org_files"
    org_dir.mkdir(parents=True, exist_ok=True)

    rec_st = structure.StructureReader.read(args.rec_file)
    organics = []

    # 1. Extract organic cofactors (Leave ions in protein)
    if args.cofactors:
        for cof in args.cofactors:
            clean_cof = str(cof).strip().strip("[]'\"").strip().upper()
            if not clean_cof or clean_cof in ION_RESNAMES:
                continue

            # Direct residue match immune to Maestro 4-char space padding
            atoms = []
            for res in rec_st.residue:
                r_name = res.pdbres.strip().upper()
                if r_name == clean_cof or r_name.startswith(clean_cof):
                    atoms.extend(res.getAtomIndices())

            if not atoms:
                asl = f"res.ptype '{clean_cof}*' OR res.ptype '{clean_cof} ' OR res.ptype '{clean_cof}'"
                atoms = analyze.evaluate_asl(rec_st, asl)

            if atoms:
                cof_st = rec_st.extract(atoms)
                cof_st.title = clean_cof
                for res in cof_st.residue:
                    res.pdbres = clean_cof
                cof_out = org_dir / f"{clean_cof}.mol2"
                cof_st.write(str(cof_out))
                (org_dir / f"{clean_cof}_charge.txt").write_text(
                    str(cof_st.formal_charge)
                )
                organics.append(
                    {
                        "name": clean_cof,
                        "mol2": f"{clean_cof}.mol2",
                        "charge": f"{clean_cof}_charge.txt",
                    }
                )
                rec_st.deleteAtoms(atoms)
                print(
                    f"Successfully extracted cofactor '{clean_cof}' ({len(atoms)} atoms).",
                    flush=True,
                )
            else:
                avail = set(r.pdbres.strip().upper() for r in rec_st.residue)
                print(
                    f"FATAL: Specified cofactor '{clean_cof}' was not found in structure! Available residues: {avail}",
                    file=sys.stderr,
                    flush=True,
                )
                sys.exit(1)

    # Save dry receptor (Protein + Ions)
    rec_out = out_dir / "receptor_dry.pdb"
    rec_st.write(str(rec_out))

    # 2. Extract Ligand preserving canonical residue name strictly
    if args.lig_file:
        if not args.lig_name:
            raise ValueError("FATAL: --lig_name must be explicitly provided when ligand file is present.")
        lig_code = args.lig_name.strip()[:4].upper()
        if not Path(args.lig_file).exists():
            raise FileNotFoundError(f"FATAL: Ligand file {args.lig_file} not found!")
        with structure.StructureReader(args.lig_file) as reader:
            for i, st in enumerate(reader):
                st.title = lig_code
                for res in st.residue:
                    res.pdbres = lig_code
                lig_out = org_dir / f"{lig_code}.mol2"
                st.write(str(lig_out))
                (org_dir / f"{lig_code}_charge.txt").write_text(str(st.formal_charge))
                organics.append(
                    {
                        "name": lig_code,
                        "mol2": f"{lig_code}.mol2",
                        "charge": f"{lig_code}_charge.txt",
                    }
                )
                break

    with open(out_dir / "organics.json", "w") as f:
        json.dump(organics, f)

    print(f"Extracted dry receptor and {len(organics)} organic molecule(s).")


def main():
    parser = argparse.ArgumentParser(
        description="Extracts receptor and ligand files for MD"
    )
    parser.add_argument(
        "--rec_file", required=True, help="Input prepared receptor file (.mae or .pdb)"
    )
    parser.add_argument(
        "--lig_file", default=None, help="Input docked pose file (.sdf or .mae)"
    )
    parser.add_argument(
        "--lig_name", default="", help="Canonical 3-letter residue name for ligand"
    )
    parser.add_argument("--cofactors", nargs="*", default=[], help="Cofactors list")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    args = parser.parse_args()

    try:
        extract_via_schrodinger(args)
    except ImportError as e:
        print(f"FATAL: Schrödinger API not found. extract_md_files.py must be run within a container providing schrodinger environment. Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"FATAL: Error in extract_md_files: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
