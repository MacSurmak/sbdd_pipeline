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
    tuple val(meta), path("index.ndx"), path("topol.top"), path("posre*.itp"), emit: restraints

    script:
    """
    echo "q" | gmx make_ndx -f ${input_gro} -o index.ndx

    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action add_groups

    # Generate multi-level restraints for Protein and Ligand (1000, 250, 50 kJ/mol/nm^2)
    echo "Protein" | gmx genrestr -f ${input_gro} -n index.ndx -o posre.itp -fc 1000 1000 1000
    echo "UNK"     | gmx genrestr -f ${input_gro} -n index.ndx -o posre_lig.itp -fc 1000 1000 1000

    echo "Protein" | gmx genrestr -f ${input_gro} -n index.ndx -o posre_250.itp -fc 250 250 250
    echo "UNK"     | gmx genrestr -f ${input_gro} -n index.ndx -o posre_lig_250.itp -fc 250 250 250

    echo "Protein" | gmx genrestr -f ${input_gro} -n index.ndx -o posre_50.itp -fc 50 50 50
    echo "UNK"     | gmx genrestr -f ${input_gro} -n index.ndx -o posre_lig_50.itp -fc 50 50 50

    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action fix_posre
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --top ${topol_top} --action inject_posre
    """

    stub:
    """
    touch index.ndx
    touch topol.top
    touch posre.itp
    touch posre_250.itp
    touch posre_50.itp
    """
}
