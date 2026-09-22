/*
 * modules/schrodinger/align_ligand.nf
 * Aligns native ligand coordinates to the PPM 3.0 oriented receptor frame (MAE & SDF)
 * Self-contained worker process: zero external script file staging.
 */

process ALIGN_LIGAND {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), path(raw_protein), path(oriented_protein), path(raw_ligand)

    output:
    tuple val(meta), path("${meta.id}_ref_native.mae"), emit: ref_native_mae
    tuple val(meta), path("${meta.id}_ref_native.sdf"), emit: ref_native_sdf

    script:
    """
    cat << 'EOF' > align_ligand.py
import argparse, sys, logging
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze, rmsd, transform

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stdout)

parser = argparse.ArgumentParser()
parser.add_argument("--raw_protein", required=True)
parser.add_argument("--oriented_protein", required=True)
parser.add_argument("--raw_ligand", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

raw_lig_path = Path(args.raw_ligand)
sdf_out = str(Path(args.output).with_suffix(".sdf"))

if not raw_lig_path.exists() or raw_lig_path.stat().st_size == 0:
    Path(args.output).touch()
    Path(sdf_out).touch()
    sys.exit(0)

raw_prot = structure.StructureReader.read(args.raw_protein, index=1)
orient_prot = structure.StructureReader.read(args.oriented_protein, index=1)
lig_st = structure.StructureReader.read(str(raw_lig_path), index=1)

calpha_asl = "atom.ptype ' CA '"
at_list_orig = analyze.evaluate_asl(raw_prot, calpha_asl)
at_list_orient = analyze.evaluate_asl(orient_prot, calpha_asl)

if len(at_list_orig) == len(at_list_orient) and len(at_list_orig) >= 3:
    matrix = rmsd.get_super_transformation_matrix(orient_prot, at_list_orient, raw_prot, at_list_orig)
    transform.transform_structure(lig_st, matrix)
    lig_st.write(args.output)
    lig_st.write(sdf_out)
else:
    lig_st.write(args.output)
    lig_st.write(sdf_out)
EOF

    schrodinger run python3 align_ligand.py \\
        --raw_protein ${raw_protein} \\
        --oriented_protein ${oriented_protein} \\
        --raw_ligand ${raw_ligand} \\
        --output ${meta.id}_ref_native.mae
    """

    stub:
    """
    touch ${meta.id}_ref_native.mae
    touch ${meta.id}_ref_native.sdf
    """
}
