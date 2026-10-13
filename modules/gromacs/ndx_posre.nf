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
    def prod_fc = (params.containsKey('prod_restraints_fc') && params.prod_restraints_fc) ? params.prod_restraints_fc : 0
    def prod_grp = (params.containsKey('prod_restraints_group') && params.prod_restraints_group) ? params.prod_restraints_group : 'C-alpha'
    def prod_fc_arg = (prod_fc.toString().toInteger() > 0) ? "--prod_fc ${prod_fc}" : ""
    def prod_grp_arg = (prod_fc.toString().toInteger() > 0) ? "--prod_group \"${prod_grp}\"" : ""
    """
    python3 ${projectDir}/bin/wrap_gromacs.py ndx_posre \\
        --input_gro ${input_gro} \\
        --topol ${topol_top} \\
        --forces "${restraints_fc}" \\
        --water_names "${water_names}" \\
        --ion_names "${ion_names}" \\
        --lipid_names "${lipid_names}" \\
        ${prod_fc_arg} \\
        ${prod_grp_arg}
    """

    stub:
    """
    touch index.ndx
    touch topol.top
    touch posre_system1.itp
    """
}