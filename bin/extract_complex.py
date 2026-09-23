#!/usr/bin/env python3
"""
bin/extract_complex.py
Separated Schrödinger worker:
Reads a macromolecular complex (.cif or .pdb), extracts protein PDB and native ligand MAE & SDF.
"""
import argparse
import sys
import logging
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze

logging.basicConfig(
    level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stdout
)


def main():
    parser = argparse.ArgumentParser(
        description="Extract components from macromolecular complex model"
    )
    parser.add_argument(
        "--input", required=True, help="Input complex file (.cif or .pdb)"
    )
    parser.add_argument(
        "--frame_index", type=int, default=1, help="NMR model index (1-based)"
    )
    parser.add_argument(
        "--ligand", required=True, help="Native ligand residue name (e.g. PKA)"
    )
    parser.add_argument(
        "--cofactors", nargs="*", default=[], help="Cofactors to retain with protein"
    )
    parser.add_argument(
        "--out_pure_protein", help="Output pure protein PDB path (for immers)"
    )
    parser.add_argument(
        "--out_complex_mae", help="Output full complex MAE path (for orientation)"
    )
    parser.add_argument(
        "--out_receptor_mae", help="Output receptor + cofactors MAE path (non-membrane)"
    )
    parser.add_argument(
        "--out_ligand_mae", help="Output native ligand MAE path (non-membrane)"
    )
    args = parser.parse_args()

    input_file = Path(args.input)
    idx = args.frame_index
    logging.info(f"Reading model {idx} from {input_file}...")

    try:
        st = structure.StructureReader.read(str(input_file), index=idx)
    except Exception as e:
        logging.error(f"Failed to read model index {idx} from {input_file}: {e}")
        sys.exit(1)

    # 1. Pure protein PDB for immers (strictly amino acids, no cofactors)
    if args.out_pure_protein:
        prot_atoms = analyze.evaluate_asl(st, "protein")
        if not prot_atoms:
            logging.error(f"No protein atoms found for model {idx} in {input_file}!")
            sys.exit(1)
        prot_st = st.extract(prot_atoms)
        prot_st.write(args.out_pure_protein)
        logging.info(
            f"Wrote pure protein PDB to {args.out_pure_protein} ({len(prot_atoms)} atoms)"
        )

    # 2. Full complex MAE for 3D transformation
    if args.out_complex_mae:
        st.write(args.out_complex_mae)
        logging.info(
            f"Wrote full complex MAE to {args.out_complex_mae} ({st.atom_total} atoms)"
        )

    # 3. Direct extraction (when membrane orientation is disabled)
    if args.out_receptor_mae:
        cofactors_asl = (
            " OR ".join([f"res.ptype '{c}'" for c in args.cofactors])
            if args.cofactors
            else ""
        )
        keep_asl = f"(protein) OR ({cofactors_asl})" if cofactors_asl else "(protein)"
        rec_atoms = analyze.evaluate_asl(st, keep_asl)
        if not rec_atoms:
            logging.error(f"No protein atoms found for model {idx} in {input_file}!")
            sys.exit(1)
        rec_st = st.extract(rec_atoms)
        rec_st.write(args.out_receptor_mae)
        logging.info(
            f"Wrote raw receptor MAE to {args.out_receptor_mae} ({len(rec_atoms)} atoms)"
        )

    if args.out_ligand_mae:
        lig_asl = f"(res.name '{args.ligand}' OR res.ptype '{args.ligand}*' OR res.ptype '*{args.ligand}*')"
        lig_atoms = analyze.evaluate_asl(st, lig_asl)
        sdf_out = str(Path(args.out_ligand_mae).with_suffix(".sdf"))
        if not lig_atoms:
            logging.warning(f"Native ligand '{args.ligand}' not found in model {idx}.")
            Path(args.out_ligand_mae).touch()
            Path(sdf_out).touch()
        else:
            lig_st = st.extract(lig_atoms)
            lig_st.write(args.out_ligand_mae)
            lig_st.write(sdf_out)
            logging.info(
                f"Wrote native ligand to {args.out_ligand_mae} and {sdf_out} ({len(lig_atoms)} atoms)"
            )


if __name__ == "__main__":
    main()
