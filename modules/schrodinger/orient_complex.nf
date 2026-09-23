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
    val(native_ligand)
    val(cofactors)
    path script_file

    output:
    tuple val(meta), path("${meta.id}_receptor_raw.mae"), emit: receptor_raw_mae
    tuple val(meta), path("${meta.id}_ref_native.mae"), emit: ref_native_mae
    tuple val(meta), path("${meta.id}_ref_native.sdf"), emit: ref_native_sdf

    script:
    def cof_arg = cofactors ? "--cofactors ${cofactors}" : ""
    """
    schrodinger run ${script_file} \\
        --raw_protein ${raw_protein} \\
        --oriented_protein ${oriented_protein} \\
        --complex_mae ${complex_mae} \\
        --ligand ${native_ligand} \\
        ${cof_arg} \\
        --out_receptor ${meta.id}_receptor_raw.mae \\
        --out_ligand ${meta.id}_ref_native.mae
    """

    stub:
    """
    touch ${meta.id}_receptor_raw.mae
    touch ${meta.id}_ref_native.mae
    touch ${meta.id}_ref_native.sdf
    """
}
