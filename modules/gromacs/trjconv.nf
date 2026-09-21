/*
 * modules/gromacs/trjconv.nf
 * PBC correction, rotational/translational fitting, and water/ion stripping
 */

process TRJCONV {
    tag "${meta.id}_rep${meta.rep}"
    label 'process_low'

    input:
    tuple val(meta), path(concat_xtc), path(tpr), path(index_ndx)

    output:
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.pdb"), emit: clean_traj

    script:
    """
    # 1. Export reference PDB of solute complex (Protein + Ligand)
    echo "PROT_LIG" | gmx trjconv -s ${tpr} -f ${tpr} -n ${index_ndx} -o complex_ref.pdb -dump 0

    # 2. PBC correction steps: whole -> center -> fit
    echo "0" | gmx trjconv -s ${tpr} -f ${concat_xtc} -o pbc_whole.xtc -pbc whole

    echo -e "PROT_LIG\\nSystem" | gmx trjconv -s ${tpr} -f pbc_whole.xtc -n ${index_ndx} -o pbc_centered.xtc -center -pbc mol -ur compact

    echo -e "Protein\\nPROT_LIG" | gmx trjconv -s ${tpr} -f pbc_centered.xtc -n ${index_ndx} -o trajectory_clean.xtc -fit rot+trans
    """

    stub:
    """
    touch trajectory_clean.xtc
    touch complex_ref.pdb
    """
}

