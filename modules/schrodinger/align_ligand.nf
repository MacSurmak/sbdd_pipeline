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
    schrodinger run python3 ${projectDir}/bin/align_ligand.py \\
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
