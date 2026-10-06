/*
 * modules/gromacs/ndx_posre.nf
 * Generates GROMACS index.ndx, thermostat groups PROT_LIG/SOL_ION, and multi-level position restraints
 */

process NDX_POSRE {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(input_gro), path(topol_top)

    output:
    tuple val(meta), path("index.ndx"), path("topol.top"), path("posre_*.itp"), emit: restraints

    script:
    """
    python3 ${projectDir}/bin/wrap_gromacs.py ndx_posre \\
        --input_gro ${input_gro} \\
        --topol ${topol_top}
    """

    stub:
    """
    touch index.ndx
    touch topol.top
    touch posre_system1.itp
    """
}