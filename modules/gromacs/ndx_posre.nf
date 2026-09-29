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
    echo "q" | gmx make_ndx -f ${input_gro} -o index.ndx
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action add_groups

    # Generate native position restraints directly from topology
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --top ${topol_top} --action generate_posre
    """

    stub:
    """
    touch index.ndx
    touch topol.top
    touch posre_system1.itp
    """
}