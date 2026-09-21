#!/usr/bin/env python3
"""
bin/align_ligand.py
Separated Schrödinger worker:
Superimposes CA atoms of raw protein to PPM-oriented protein,
applies the transformation matrix to native ligand, and saves aligned ligand (MAE & SDF).
"""
import argparse
import sys
import logging
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze, rmsd, transform

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stdout)

def main():
    parser = argparse.ArgumentParser(description="Align native ligand to PPM-oriented receptor frame")
    parser.add_argument("--raw_protein", required=True, help="Input un-oriented protein PDB")
    parser.add_argument("--oriented_protein", required=True, help="Input PPM-oriented protein PDB")
    parser.add_argument("--raw_ligand", required=True, help="Input raw native ligand MAE")
    parser.add_argument("--output", required=True, help="Output aligned native ligand MAE")
    args = parser.parse_args()

    raw_lig_path = Path(args.raw_ligand)
    sdf_out = str(Path(args.output).with_suffix(".sdf"))

    if not raw_lig_path.exists() or raw_lig_path.stat().st_size == 0:
        logging.warning("No valid raw ligand file found. Writing empty placeholders.")
        Path(args.output).touch()
        Path(sdf_out).touch()
        sys.exit(0)

    try:
        raw_prot = structure.StructureReader.read(args.raw_protein, index=1)
        orient_prot = structure.StructureReader.read(args.oriented_protein, index=1)
        lig_st = structure.StructureReader.read(str(raw_lig_path), index=1)
    except Exception as e:
        logging.error(f"Failed to read input structures for alignment: {e}")
        sys.exit(1)

    calpha_asl = "atom.ptype ' CA '"
    at_list_orig = analyze.evaluate_asl(raw_prot, calpha_asl)
    at_list_orient = analyze.evaluate_asl(orient_prot, calpha_asl)

    if len(at_list_orig) == len(at_list_orient) and len(at_list_orig) >= 3:
        matrix = rmsd.get_super_transformation_matrix(orient_prot, at_list_orient, raw_prot, at_list_orig)
        transform.transform_structure(lig_st, matrix)
        lig_st.write(args.output)
        lig_st.write(sdf_out)
        logging.info(f"Successfully aligned native ligand -> {args.output} and {sdf_out}")
    else:
        logging.warning("CA atom mismatch during alignment; writing untransformed ligand.")
        lig_st.write(args.output)
        lig_st.write(sdf_out)

if __name__ == "__main__":
    main()
