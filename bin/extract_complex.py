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

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stdout)

def main():
    parser = argparse.ArgumentParser(description="Extract protein PDB and native ligand from complex")
    parser.add_argument("--input", required=True, help="Input complex file (.cif or .pdb)")
    parser.add_argument("--frame_index", type=int, default=1, help="NMR model index (1-based)")
    parser.add_argument("--ligand", required=True, help="Native ligand residue name (e.g. PKA)")
    parser.add_argument("--cofactors", nargs="*", default=[], help="Cofactors to retain with protein")
    parser.add_argument("--out_protein", required=True, help="Output protein PDB path")
    parser.add_argument("--out_ligand", required=True, help="Output native ligand MAE path")
    args = parser.parse_args()

    input_file = Path(args.input)
    idx = args.frame_index
    logging.info(f"Reading model {idx} from {input_file}...")

    try:
        st = structure.StructureReader.read(str(input_file), index=idx)
    except Exception as e:
        logging.error(f"Failed to read model index {idx} from {input_file}: {e}")
        sys.exit(1)

    # 1. Extract protein & cofactors
    cofactors_asl = " OR ".join([f"res.ptype '{c}'" for c in args.cofactors])
    keep_asl = f"(protein) OR ({cofactors_asl})" if cofactors_asl else "(protein)"
    rec_atoms = analyze.evaluate_asl(st, keep_asl)
    if not rec_atoms:
        logging.error(f"No protein atoms found for model {idx} in {input_file}!")
        sys.exit(1)

    rec_st = st.extract(rec_atoms)
    rec_st.write(args.out_protein)
    logging.info(f"Wrote raw protein to {args.out_protein} ({len(rec_atoms)} atoms)")

    # 2. Extract native ligand (MAE and SDF)
    lig_asl = f"res.ptype '{args.ligand}'"
    lig_atoms = analyze.evaluate_asl(st, lig_asl)
    sdf_out = str(Path(args.out_ligand).with_suffix(".sdf"))
    if not lig_atoms:
        logging.warning(f"Native ligand '{args.ligand}' not found in model {idx}.")
        Path(args.out_ligand).touch()
        Path(sdf_out).touch()
    else:
        lig_st = st.extract(lig_atoms)
        lig_st.write(args.out_ligand)
        lig_st.write(sdf_out)
        logging.info(f"Wrote native ligand to {args.out_ligand} and {sdf_out} ({len(lig_atoms)} atoms)")

if __name__ == "__main__":
    main()
