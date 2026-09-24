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
    parser.add_argument("--ligand", required=True, help="Native ligand residue name")
    parser.add_argument(
        "--cofactors", nargs="*", default=[], help="Cofactors to retain with receptor"
    )
    parser.add_argument(
        "--out_receptor", required=True, help="Output oriented receptor+cofactors MAE"
    )
    parser.add_argument(
        "--out_ligand", required=True, help="Output oriented native ligand MAE"
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

    if len(at_list_orig) == len(at_list_orient) and len(at_list_orig) >= 3:
        matrix = rmsd.get_super_transformation_matrix(
            orient_prot, at_list_orient, raw_prot, at_list_orig
        )
        transform.transform_structure(complex_st, matrix)
        logging.info("Successfully applied PPM 3.0 transformation to the full complex.")
    else:
        logging.warning(
            "CA atom mismatch during alignment; using untransformed complex coordinates."
        )

    # 1. Extract native ligand (MAE & SDF)
    lig_asl = f"res.ptype '{args.ligand}'"
    lig_atoms = analyze.evaluate_asl(complex_st, lig_asl)
    sdf_out = str(Path(args.out_ligand).with_suffix(".sdf"))
    if not lig_atoms:
        logging.error(
            f"STRICT MODE: Native ligand '{args.ligand}' not found in complex."
        )
        sys.exit(1)

    lig_st = complex_st.extract(lig_atoms)
    lig_st.write(args.out_ligand)
    lig_st.write(sdf_out)
    logging.info(f"Wrote oriented native ligand to {args.out_ligand} and {sdf_out}")
    # 2. Extract receptor + cofactors (MAE)
    cofactors_asl = (
        " OR ".join([f"res.ptype '{c}'" for c in args.cofactors])
        if args.cofactors
        else ""
    )
    keep_asl = f"(protein) OR ({cofactors_asl})" if cofactors_asl else "(protein)"
    rec_atoms = analyze.evaluate_asl(complex_st, keep_asl)
    if not rec_atoms:
        logging.error("No protein atoms found in oriented complex!")
        sys.exit(1)

    rec_st = complex_st.extract(rec_atoms)
    rec_st.write(args.out_receptor)
    logging.info(
        f"Wrote oriented raw receptor to {args.out_receptor} ({len(rec_atoms)} atoms)"
    )


if __name__ == "__main__":
    main()
