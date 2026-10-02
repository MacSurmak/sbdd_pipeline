/*
 * modules/gnina/rescore.nf
 * GNINA CNN Rescoring / Minimization Module
 */

process GNINA_RESCORE {
    tag "${meta.id}"
    label 'process_gpu'

    input:
    tuple val(meta), path(receptor_pdb, stageAs: 'rec_dir/*'), path(ref_lig_sdf, stageAs: 'ref_dir/*'), path(poses_sdf, stageAs: 'poses_dir/*')
    val(cnn_model)

    output:
    tuple val(meta), path("${meta.id}_rescored.sdf"), emit: rescored_sdf
    tuple val(meta), path("${meta.id}_gnina_rescore.log"), emit: log

    script:
    def cnn_flag = cnn_model ? "--cnn ${cnn_model}" : ""
    """
    if [ ! -f "${ref_lig_sdf}" ] || [ ! -s "${ref_lig_sdf}" ]; then
        echo "FATAL: Reference ligand SDF is missing or empty: ${ref_lig_sdf}"
        exit 1
    fi
    AUTOBOX_FLAG="--autobox_ligand ${ref_lig_sdf}"

    gnina \\
        -r ${receptor_pdb} \\
        -l ${poses_sdf} \\
        \$AUTOBOX_FLAG \\
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

