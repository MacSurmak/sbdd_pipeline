#!/usr/bin/env python3
"""
bin/split_prepared_complex.py
Splits a fully prepared and protonated complex (from PrepWizard) into
receptor (protein+cofactors) and native ligand.
Preserves perfect in-pocket coordinates and correctly assigned protonation states.
"""
import argparse
import sys
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze


def main():
    parser = argparse.ArgumentParser(
        description="Split prepared complex into receptor and ligand"
    )
    parser.add_argument(
        "--complex_mae", required=True, help="Input prepared complex MAE"
    )
    parser.add_argument("--ligand", default="", help="Native ligand residue name")
    parser.add_argument(
        "--out_receptor_mae", required=True, help="Output prepared receptor MAE"
    )
    parser.add_argument(
        "--out_receptor_pdb", required=True, help="Output prepared receptor PDB"
    )
    parser.add_argument(
        "--out_ligand_mae", required=True, help="Output prepared ligand MAE"
    )
    parser.add_argument(
        "--out_ligand_sdf", required=True, help="Output prepared ligand SDF"
    )
    args = parser.parse_args()

    try:
        with structure.StructureReader(args.complex_mae) as reader:
            sts = list(reader)
        if not sts:
            raise ValueError(f"No structures found in {args.complex_mae}")
        st = sts[0]
        for next_st in sts[1:]:
            st.extend(next_st)
    except Exception as e:
        print(f"FATAL: Could not read {args.complex_mae}: {e}")
        sys.exit(1)

    if args.ligand:
        target_res = args.ligand.strip().upper()
        lig_atoms = []

        # 1. Direct Python residue match immune to Maestro 4-char space padding
        for res in st.residue:
            r_name = res.pdbres.strip().upper()
            if r_name == target_res or r_name.startswith(target_res):
                lig_atoms.extend(res.getAtomIndices())

        # 2. Fallback ASL with wildcards
        if not lig_atoms:
            lig_atoms = analyze.evaluate_asl(
                st, f"res.ptype '{target_res}*' OR res.ptype '{target_res} ' OR res.ptype '{target_res}'"
            )

        if lig_atoms:
            lig_st = st.extract(lig_atoms)
            lig_st.title = args.ligand
            lig_st.write(args.out_ligand_mae)
            lig_st.write(args.out_ligand_sdf)
            st.deleteAtoms(lig_atoms)
            print(
                f"Extracted prepared native ligand '{args.ligand}' ({len(lig_atoms)} atoms)."
            )
        else:
            avail = set(r.pdbres.strip() for r in st.residue)
            raise ValueError(
                f"FATAL: Native ligand '{args.ligand}' not found in prepared complex! "
                f"Available residues: {avail}"
            )
    else:
        # If no ligand was specified, write pure receptor and exit without creating dummy 0-byte ligand files
        st.write(args.out_receptor_mae)
        st.write(args.out_receptor_pdb)
        print(f"Wrote prepared apo receptor to {args.out_receptor_mae} ({st.atom_total} atoms).")
        return

    # The remaining structure is Protein + Cofactors
    st.write(args.out_receptor_mae)
    st.write(args.out_receptor_pdb)
    print(
        f"Wrote prepared receptor to {args.out_receptor_mae} and {args.out_receptor_pdb} ({st.atom_total} atoms)."
    )


if __name__ == "__main__":
    main()
