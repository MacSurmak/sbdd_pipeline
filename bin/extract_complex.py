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
        "--ligand", default="", help="Native ligand residue name (e.g. PKA)"
    )
    parser.add_argument(
        "--cofactors", nargs="*", default=[], help="Cofactors to retain with protein"
    )
    parser.add_argument(
        "--out_pure_protein", help="Output pure protein PDB path (for immers)"
    )
    parser.add_argument(
        "--out_complex_mae", required=True, help="Output full complex MAE path"
    )
    args = parser.parse_args()

    input_file = Path(args.input)
    idx = args.frame_index
    logging.info(f"Reading model {idx} from {input_file}...")

    try:
        with structure.StructureReader(str(input_file)) as reader:
            sts = list(reader)

        # If it is a multi-CT complex (Receptor CT + Ligand CT), merge all CTs into one unified structure
        if len(sts) > 1 and any(analyze.evaluate_asl(s, "protein") for s in sts) and any(analyze.evaluate_asl(s, "not protein and not water") for s in sts) and not all(s.atom_total == sts[0].atom_total for s in sts):
            logging.info("Detected multi-CT complex in MAE file. Merging all CT blocks into unified complex.")
            complex_st = sts[0]
            for other_st in sts[1:]:
                complex_st.extend(other_st)
            st = complex_st
        else:
            st = sts[idx - 1]
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

    # 2. Write Full Complex MAE
    st.write(args.out_complex_mae)
    logging.info(
        f"Wrote full complex MAE to {args.out_complex_mae} ({st.atom_total} atoms)"
    )


if __name__ == "__main__":
    main()
