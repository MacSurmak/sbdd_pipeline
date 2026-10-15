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
    tuple val(meta), path("rescored.sdf"), emit: rescored_sdf
    tuple val(meta), path("rescored.log"), emit: log

    script:
    def cnn_arg = cnn_model ? "--cnn_model ${cnn_model}" : ""
    """
    trap '[[ -f rescored.log ]] && cp -f rescored.log ${task.workDir}/ 2>/dev/null || true' EXIT ERR
    python3 ${projectDir}/bin/wrap_gnina.py rescore \\
        --receptor ${receptor_pdb} \\
        --ligand ${poses_sdf} \\
        --ref_ligand ${ref_lig_sdf} \\
        ${cnn_arg} \\
        --out_sdf rescored.sdf \\
        --out_log rescored.log
    """
}

