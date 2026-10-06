/*
 * modules/gromacs/ndx_posre.nf
 * Generates GROMACS index.ndx, thermostat groups PROT_LIG/SOL_ION, and multi-level position restraints
 */

process NDX_POSRE {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(input_gro), path(topol_top)
    val(restraints_fc)
    val(water_names)
    val(ion_names)
    val(lipid_names)

    output:
    tuple val(meta), path("index.ndx"), path("topol.top"), path("posre_*.itp"), emit: restraints

    script:
    """
    python3 ${projectDir}/bin/wrap_gromacs.py ndx_posre \\
        --input_gro ${input_gro} \\
        --topol ${topol_top} \\
        --forces "${restraints_fc}" \\
        --water_names "${water_names}" \\
        --ion_names "${ion_names}" \\
        --lipid_names "${lipid_names}"
    """

    stub:
    """
    touch index.ndx
    touch topol.top
    touch posre_system1.itp
    """
}