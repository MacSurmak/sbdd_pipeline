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
    path mdp_templates
    val(is_membrane)
    val(temperature)
    val(total_ns)
    val(chunk_ns)
    val(slurm_queue)

    output:
    tuple val(meta), path("md_*.xtc"), emit: all_xtcs
    tuple val(meta), path("md_1.tpr"), path(topol_top), path(index_ndx), emit: sys_files

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
        --mdp_template md.mdp \\
        --cpus ${task.cpus} \\
        --outdir .
    """

    stub:
    """
    touch md_1.xtc
    touch md_1.tpr
    """
}
