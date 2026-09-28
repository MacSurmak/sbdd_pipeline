/*
 * modules/gromacs/trjconv.nf
 * PBC correction, rotational/translational fitting, and water/ion stripping
 */

process TRJCONV {
    tag "${meta.id}_rep${meta.rep}"
    label 'process_low'

    input:
    tuple val(meta), path(concat_xtc), path(tpr), path(topol_top), path(index_ndx)
    val(pbc_center)
    val(pbc_extract)

    output:
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.pdb"), path(topol_top), path(index_ndx), emit: clean_traj

    script:
    """
    # 1. Inject [ COMPLEX ] (Protein + UNK strictly, no lipids) into index.ndx if not present
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx ${index_ndx} --action ensure_complex

    # 2. PBC correction steps: whole -> center -> extract selected group
    echo "0" | gmx trjconv -s ${tpr} -f ${concat_xtc} -o pbc_whole.xtc -pbc whole

    echo -e "${pbc_center}\\n${pbc_extract}" | gmx trjconv -s ${tpr} -f pbc_whole.xtc -n ${index_ndx} -o trajectory_clean.xtc -center -pbc mol -ur compact

    echo "${pbc_extract}" | gmx trjconv -s ${tpr} -f trajectory_clean.xtc -n ${index_ndx} -o complex_ref.pdb -dump 0

    rm -f pbc_whole.xtc
    """

    stub:
    """
    touch trajectory_clean.xtc
    touch complex_ref.pdb
    """
}

