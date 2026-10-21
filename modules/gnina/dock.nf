/*
 * modules/gnina/dock.nf
 * Distributed GNINA Docking Process for Slurm GPU nodes
 */

process GNINA_DOCK {
    tag "${meta.id}_chunk_${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), path(receptor_pdb, stageAs: 'rec_dir/*'), path(ref_lig_sdf, stageAs: 'ref_dir/*'), path(chunk_sdf, stageAs: 'chunk_dir/*'), val(chunk_idx)
    val(exhaustiveness)
    val(cnn_scoring)
    val(num_modes)
    val(autobox_add)
    output:
    tuple val(meta), path("gnina_dock.sdf"), emit: docked_sdf
    tuple val(meta), path("gnina_dock.log"), emit: log

    script:
    """
    python3 ${projectDir}/bin/wrap_gnina.py dock \\
        --receptor ${receptor_pdb} \\
        --ligand ${chunk_sdf} \\
        --ref_ligand ${ref_lig_sdf} \\
        --autobox_add ${autobox_add} \\
        --exhaustiveness ${exhaustiveness} \\
        --cnn_scoring ${cnn_scoring} \\
        --num_modes ${num_modes} \\
        --out_sdf gnina_dock.sdf \\
        --out_log gnina_dock.log
    """
}

