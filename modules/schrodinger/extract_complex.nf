/*
 * modules/schrodinger/extract_complex.nf
 * Extracts raw protein PDB and native ligand MAE & SDF from complex (.cif / .pdb)
 * Self-contained worker process: zero external script file staging.
 */

process EXTRACT_COMPLEX {
    tag "${meta.id}_f${frame_idx}"
    label "process_single"

    input:
    tuple val(meta), path(complex_file), val(frame_idx)
    val(native_ligand)
    val(cofactors)

    output:
    tuple val(meta), path("${meta.id}_raw_protein.pdb"), emit: raw_protein_pdb
    tuple val(meta), path("${meta.id}_raw_native.mae"), emit: raw_native_mae
    tuple val(meta), path("${meta.id}_raw_native.sdf"), emit: raw_native_sdf

    script:
    def cof_arg = cofactors ? "--cofactors ${cofactors}" : ""
    """
    cat << 'EOF' > extract_complex.py
import argparse, sys, logging
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stdout)

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--frame_index", type=int, default=1)
parser.add_argument("--ligand", required=True)
parser.add_argument("--cofactors", nargs="*", default=[])
parser.add_argument("--out_protein", required=True)
parser.add_argument("--out_ligand", required=True)
args = parser.parse_args()

st = structure.StructureReader.read(args.input, index=args.frame_index)
cofactors_asl = " OR ".join([f"res.ptype '{c}'" for c in args.cofactors])
keep_asl = f"(protein) OR ({cofactors_asl})" if cofactors_asl else "(protein)"
rec_atoms = analyze.evaluate_asl(st, keep_asl)
if not rec_atoms:
    logging.error("No protein atoms found")
    sys.exit(1)

rec_st = st.extract(rec_atoms)
rec_st.write(args.out_protein)

lig_asl = f"res.ptype '{args.ligand}'"
lig_atoms = analyze.evaluate_asl(st, lig_asl)
sdf_out = str(Path(args.out_ligand).with_suffix(".sdf"))
if not lig_atoms:
    Path(args.out_ligand).touch()
    Path(sdf_out).touch()
else:
    lig_st = st.extract(lig_atoms)
    lig_st.write(args.out_ligand)
    lig_st.write(sdf_out)
EOF

    schrodinger run python3 extract_complex.py \\
        --input ${complex_file} \\
        --frame_index ${frame_idx} \\
        --ligand ${native_ligand} \\
        ${cof_arg} \\
        --out_protein ${meta.id}_raw_protein.pdb \\
        --out_ligand ${meta.id}_raw_native.mae
    """

    stub:
    """
    touch ${meta.id}_raw_protein.pdb
    touch ${meta.id}_raw_native.mae
    touch ${meta.id}_raw_native.sdf
    """
}
