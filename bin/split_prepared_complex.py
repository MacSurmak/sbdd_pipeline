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
        st = structure.StructureReader.read(args.complex_mae)
    except Exception as e:
        print(f"FATAL: Could not read {args.complex_mae}: {e}")
        sys.exit(1)

    if args.ligand:
        lig_atoms = analyze.evaluate_asl(st, f"res.ptype '{args.ligand}'")
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
            print(
                f"Warning: Native ligand '{args.ligand}' not found in prepared complex."
            )
            Path(args.out_ligand_mae).touch()
            Path(args.out_ligand_sdf).touch()
    else:
        Path(args.out_ligand_mae).touch()
        Path(args.out_ligand_sdf).touch()

    # The remaining structure is Protein + Cofactors
    st.write(args.out_receptor_mae)
    st.write(args.out_receptor_pdb)
    print(
        f"Wrote prepared receptor to {args.out_receptor_mae} and {args.out_receptor_pdb} ({st.atom_total} atoms)."
    )


if __name__ == "__main__":
    main()
