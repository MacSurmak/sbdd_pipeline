/*
 * modules/amber/parmed.nf
 * Converts Amber prmtop/inpcrd to GROMACS topol.top and step3_input.gro
 */

process PARMED {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(prmtop), path(inpcrd)
    val(hmr_enabled)
    val(hmr_mass)

    output:
    tuple val(meta), path("topol.top"), path("step3_input.gro"), emit: gmx_system

    script:
    def hmr_arg = (hmr_enabled.toString() == 'true') ? "--hmr --hmr_mass ${hmr_mass}" : ""
    """
    python3 ${projectDir}/bin/convert_amber_to_gmx.py --prmtop ${prmtop} --inpcrd ${inpcrd} ${hmr_arg}
    """
}

