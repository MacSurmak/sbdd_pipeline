/*
 * modules/unidock/dock.nf
 * Distributed GPU-Accelerated Uni-Dock Sampling Process for Slurm GPU nodes
 */

process UNIDOCK_DOCK {
    tag "${meta.id}"
    label 'process_unidock'

    input:
    tuple val(meta), path(receptor_pdbqt, stageAs: 'rec_dir/*'), path(ref_lig_sdf, stageAs: 'ref_dir/*'), path(chunk_sdf, stageAs: 'chunk_dir/*'), val(chunk_idx)
    val(search_mode)
    val(scoring)
    val(refine_step)
    val(num_modes)
    val(energy_range)
    val(autobox_add)
    val(seed)

    output:
    tuple val(meta), path("unidock_dock.sdf"), emit: docked_sdf
    tuple val(meta), path("unidock_dock.log"), emit: log

    script:
    """
    python3 ${projectDir}/bin/wrap_unidock.py \\
        --receptor ${receptor_pdbqt} \\
        --ligand ${chunk_sdf} \\
        --ref_ligand ${ref_lig_sdf} \\
        --search_mode ${search_mode} \\
        --scoring ${scoring} \\
        --refine_step ${refine_step} \\
        --num_modes ${num_modes} \\
        --energy_range ${energy_range} \\
        --autobox_add ${autobox_add} \\
        --seed ${seed} \\
        --out_sdf unidock_dock.sdf \\
        --out_log unidock_dock.log
    """
}
