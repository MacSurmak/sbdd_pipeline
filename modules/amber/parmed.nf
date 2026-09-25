/*
 * modules/amber/parmed.nf
 * Converts Amber prmtop/inpcrd to GROMACS topol.top and step3_input.gro
 */

process PARMED {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(prmtop), path(inpcrd)

    output:
    tuple val(meta), path("topol.top"), path("step3_input.gro"), emit: gmx_system

    script:
    """
    convert_amber_to_gmx.py --prmtop ${prmtop} --inpcrd ${inpcrd}
    """

    stub:
    """
    touch topol.top
    touch step3_input.gro
    """
}

