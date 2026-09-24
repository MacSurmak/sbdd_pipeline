/*
 * modules/gnina/rescore.nf
 * GNINA CNN Rescoring / Minimization Module
 */

process GNINA_RESCORE {
    tag "${meta.id}"
    label 'process_gpu'

    input:
    tuple val(meta), path(receptor_pdb), path(ref_lig_sdf), path(poses_sdf)
    val(cnn_model)

    output:
    tuple val(meta), path("${meta.id}_rescored.sdf"), emit: rescored_sdf
    tuple val(meta), path("${meta.id}_gnina_rescore.log"), emit: log

    script:
    def cnn_flag = cnn_model ? "--cnn ${cnn_model}" : ""
    def autobox_flag = (ref_lig_sdf && ref_lig_sdf.name != 'NO_FILE') ? "--autobox_ligand ${ref_lig_sdf}" : ""
    """
    gnina \\
        -r ${receptor_pdb} \\
        -l ${poses_sdf} \\
        ${autobox_flag} \\
        --score_only \\
        ${cnn_flag} \\
        -o ${meta.id}_rescored.sdf > ${meta.id}_gnina_rescore.log
    """

    stub:
    """
    touch ${meta.id}_rescored.sdf
    touch ${meta.id}_gnina_rescore.log
    """
}

