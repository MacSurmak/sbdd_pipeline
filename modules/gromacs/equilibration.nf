/*
 * modules/gromacs/equilibration.nf
 * Executes multi-stage energy minimization and equilibration (NVT/NPT) on Slurm GPU
 */

process GMX_EQUILIBRATION {
    tag "${meta.id}"
    label 'process_gpu'

    input:
    tuple val(meta), path(input_gro), path(topol_top), path(index_ndx), path(posre_files)
    path md_protocol
    val minim_flags
    val eq_flags

    output:
    tuple val(meta), path("eq_final.gro"), path("eq_final.cpt"), path("eq_final.tpr"), emit: equilibrated_state
    tuple val(meta), path("*.log"), emit: logs

    script:
    """
    python3 ${projectDir}/bin/wrap_gromacs.py equilibration \\
        --input_gro ${input_gro} \\
        --topol ${topol_top} \\
        --ndx ${index_ndx} \\
        --cpus ${task.cpus} \\
        --posre_files posre_*.itp \\
        --sys_name "${meta.id}" \\
        --minim_flags "${minim_flags}" \\
        --eq_flags "${eq_flags}"
    """

    stub:
    """
    touch eq_final.gro
    touch eq_final.cpt
    touch eq_final.tpr
    touch eq.log
    """
}
