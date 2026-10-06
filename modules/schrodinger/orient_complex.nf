/*
 * modules/schrodinger/orient_complex.nf
 * Transforms full complex structure using PPM 3.0 CA superposition matrix.
 * Extracts oriented receptor+cofactors (MAE) and native ligand (MAE & SDF).
 */
process ORIENT_COMPLEX {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), path(raw_protein), path(oriented_protein), path(complex_mae)

    output:
    tuple val(meta), path("complex_oriented.mae"), emit: oriented_mae

    script:
    """
    schrodinger run ${projectDir}/bin/orient_complex.py \\
        --raw_protein ${raw_protein} \\
        --oriented_protein ${oriented_protein} \\
        --complex_mae ${complex_mae} \\
        --out_complex_mae complex_oriented.mae
    """
    stub:
    """
    touch ${meta.id}_complex_oriented.mae
    """
}
