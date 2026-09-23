/*
 * modules/gromacs/run_chunk.nf
 * Orchestrates dynamic atomic production MD simulation chunks (1 to N chunks)
 * Completely eliminates any hardcoded chunk limits or manual DAG unrolling.
 */

process RUN_PRODUCTION_MD {
    tag "${meta.id}_rep${meta.rep}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(eq_gro), path(eq_cpt), path(topol_top), path(index_ndx)
    val(is_membrane)
    val(temperature)
    val(total_ns)
    val(chunk_ns)
    val(slurm_queue)

    output:
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.pdb"), path(topol_top), path(index_ndx), emit: clean_traj
    tuple val(meta), path("md_*.xtc"), emit: all_xtcs

    script:
    """
    run_production_md.py \\
        --meta_id ${meta.id} \\
        --rep ${meta.rep} \\
        --gro ${eq_gro} \\
        --cpt ${eq_cpt} \\
        --top ${topol_top} \\
        --ndx ${index_ndx} \\
        --total_ns ${total_ns} \\
        --chunk_ns ${chunk_ns} \\
        --temp ${temperature} \\
        --is_membrane ${is_membrane} \\
        --queue ${slurm_queue} \\
        --outdir .
    """

    stub:
    """
    touch md_1.xtc
    touch trajectory_clean.xtc
    touch complex_ref.pdb
    """
}
