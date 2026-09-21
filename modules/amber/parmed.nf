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
    python3 -c "import parmed as p; s=p.load_file('${prmtop}', '${inpcrd}'); s.save('topol.top', format='gromacs'); s.save('step3_input.gro')"
    """

    stub:
    """
    touch topol.top
    touch step3_input.gro
    """
}

