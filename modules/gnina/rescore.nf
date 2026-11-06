/*
 * modules/gnina/rescore.nf
 * GNINA CNN Rescoring / Minimization Module
 */

process GNINA_RESCORE {
    tag "${meta.id}_${engine_tag}"
    label 'process_gpu_shard'

    input:
    tuple val(meta), path(receptor_pdb, stageAs: 'rec_dir/*'), path(ref_lig_sdf, stageAs: 'ref_dir/*'), path(poses_sdf, stageAs: 'poses_dir/*')
    val(engine_tag)
    val(rescore_mode)
    val(autobox_add)
    val(cnn_model)

    output:
    tuple val(meta), path("${meta.id}_${engine_tag}_rescored.sdf"), emit: rescored_sdf
    tuple val(meta), path("${meta.id}_${engine_tag}_rescored.log"), emit: log

    script:
    """
    python3 ${projectDir}/bin/wrap_gnina.py rescore \\
        --receptor ${receptor_pdb} \\
        --ligand ${poses_sdf} \\
        --ref_ligand ${ref_lig_sdf} \\
        --cpus ${task.cpus} \\
        --rescore_mode ${rescore_mode} \\
        --autobox_add ${autobox_add} \\
        --cnn_model "${cnn_model}" \\
        --out_sdf "${meta.id}_${engine_tag}_rescored.sdf" \\
        --out_log "${meta.id}_${engine_tag}_rescored.log"
    """
}

