#!/usr/bin/env python3
"""
bin/orient_complex.py
Calculates PPM 3.0 rotation/translation matrix between unoriented pure protein PDB
and oriented protein PDB (from immers), applies the transformation to the entire
Schrodinger complex MAE, and extracts:
1. Native ligand: ref_native.mae and ref_native.sdf
2. Receptor + cofactors: receptor_raw.mae
Preserves all chemistry, bond orders, formal charges, and cofactors.
"""
import argparse
import sys
import logging
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze, rmsd, transform

logging.basicConfig(
    level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stdout
)


def main():
    parser = argparse.ArgumentParser(
        description="Orient complex structure via PPM 3.0 CA superposition"
    )
    parser.add_argument(
        "--raw_protein", required=True, help="Input un-oriented pure protein PDB"
    )
    parser.add_argument(
        "--oriented_protein", required=True, help="Input PPM-oriented protein PDB"
    )
    parser.add_argument("--complex_mae", required=True, help="Input raw complex MAE")
    parser.add_argument(
        "--out_complex_mae", required=True, help="Output oriented full complex MAE"
    )
    args = parser.parse_args()

    try:
        raw_prot = structure.StructureReader.read(args.raw_protein, index=1)
        orient_prot = structure.StructureReader.read(args.oriented_protein, index=1)
        complex_st = structure.StructureReader.read(args.complex_mae, index=1)
    except Exception as e:
        logging.error(f"Failed to read input structures for orientation: {e}")
        sys.exit(1)

    calpha_asl = "atom.ptype ' CA '"
    at_list_orig = analyze.evaluate_asl(raw_prot, calpha_asl)
    at_list_orient = analyze.evaluate_asl(orient_prot, calpha_asl)

    if len(at_list_orig) != len(at_list_orient) or len(at_list_orig) < 3:
        logging.error(
            "FATAL: CA atom mismatch between original and PPM-oriented protein. Alignment failed."
        )
        sys.exit(1)

    matrix = rmsd.get_super_transformation_matrix(
        orient_prot, at_list_orient, raw_prot, at_list_orig
    )
    transform.transform_structure(complex_st, matrix)
    logging.info("Successfully applied PPM 3.0 transformation to the full complex.")

    complex_st.write(args.out_complex_mae)
    logging.info(f"Wrote oriented full complex to {args.out_complex_mae}")


if __name__ == "__main__":
    main()
