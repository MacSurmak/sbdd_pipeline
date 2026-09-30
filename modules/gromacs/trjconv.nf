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
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.pdb"), path(topol_top, includeInputs: true), path(index_ndx, includeInputs: true), emit: clean_traj

    script:
    """
    # 1. Build strict [ COMPLEX ] group (Protein + Ligand + Cofactors, stripped of water and ions)
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx ${index_ndx} --action ensure_complex

    # 2. PBC Step 1: Make all molecules whole across periodic boundaries
    echo "0" | gmx trjconv -s ${tpr} -f ${concat_xtc} -o pbc_whole.xtc -pbc whole

    # 3. PBC Step 2: Center on Protein, make molecules compact, and extract full COMPLEX
    echo -e "${pbc_center}\\n${pbc_extract}" | gmx trjconv -s ${tpr} -f pbc_whole.xtc -n ${index_ndx} -o pbc_centered.xtc -center -pbc mol -ur compact

    # 4. PBC Step 3: Rotational and translational fitting on Protein backbone to eliminate drift
    echo -e "Protein\\n${pbc_extract}" | gmx trjconv -s ${tpr} -f pbc_centered.xtc -n ${index_ndx} -o trajectory_clean.xtc -fit rot+trans

    # 5. Extract reference 3D structure of the dry COMPLEX at t=0
    echo "${pbc_extract}" | gmx trjconv -s ${tpr} -f trajectory_clean.xtc -n ${index_ndx} -o complex_ref.pdb -dump 0

    rm -f pbc_whole.xtc pbc_centered.xtc
    """

    stub:
    """
    touch trajectory_clean.xtc
    touch complex_ref.pdb
    """
}

