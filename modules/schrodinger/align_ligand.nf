/*
 * modules/schrodinger/align_ligand.nf
 * Aligns native ligand coordinates to the PPM 3.0 oriented receptor frame (MAE & SDF)
 */

process ALIGN_LIGAND {
    tag "${meta.id}"
    label "process_single"
    label "process_master"

    input:
    tuple val(meta), path(raw_protein), path(oriented_protein), path(raw_ligand)
    path(align_script)

    output:
    tuple val(meta), path("${meta.id}_ref_native.mae"), emit: ref_native_mae
    tuple val(meta), path("${meta.id}_ref_native.sdf"), emit: ref_native_sdf

    script:
    """
    schrodinger run python3 ${align_script} \\
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
