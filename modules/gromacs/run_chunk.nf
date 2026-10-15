/*
 * modules/gromacs/run_chunk.nf
 * Orchestrates dynamic atomic production MD simulation chunks (1 to N chunks)
 * Completely eliminates any hardcoded chunk limits or manual DAG unrolling.
 */

process RUN_MD_CHUNK {
    tag "${meta.id}_R${meta.rep}_C${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), path(prev_gro), path(prev_cpt), path(prev_tpr), path(topol_top), path(index_ndx), path(posre_files)
    path md_protocol
    val(chunk_ns)
    val(chunk_idx)
    val prod_flags

    output:
    tuple val(meta), path("out.gro"), path("out.cpt"), path("out.tpr"), path(topol_top, includeInputs: true), path(index_ndx, includeInputs: true), path("posre_*.itp", includeInputs: true), emit: state
    tuple val(meta), path("md_${chunk_idx}.xtc"), emit: xtc
    tuple val(meta), path("sys_ref.tpr"), path(topol_top, includeInputs: true), path(index_ndx, includeInputs: true), emit: sys_files

    script:
    """
    python3 ${projectDir}/bin/wrap_gromacs.py run_chunk \\
        --prev_gro ${prev_gro} \\
        --prev_cpt ${prev_cpt} \\
        --prev_tpr ${prev_tpr} \\
        --topol ${topol_top} \\
        --ndx ${index_ndx} \\
        --chunk_ns ${chunk_ns} \\
        --chunk_idx ${chunk_idx} \\
        --cpus ${task.cpus} \\
        --posre_files posre_*.itp \\
        --sys_name "${meta.id}_R${meta.rep}_C${chunk_idx}" \\
        --prod_flags "${prod_flags}"
    """
}
