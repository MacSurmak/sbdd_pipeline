/*
 * modules/gromacs/trjconv.nf
 * PBC correction, rotational/translational fitting, and water/ion stripping
 */

process TRJCONV {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(xtc_inputs, stageAs: 'inputs/*'), path(tpr)
    val(pbc_center)
    val(pbc_extract)

    output:
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.gro"), emit: clean_traj

    script:
    """
    # 1. Concatenate chunks strictly in worker scratch without duplicating on project storage
    XTC_COUNT=\$(ls -1 inputs/*.xtc 2>/dev/null | wc -l)
    if [ "\$XTC_COUNT" -gt 1 ]; then
        ls -v inputs/*.xtc > sorted_chunks.list
        gmx trjcat -f \$(cat sorted_chunks.list) -o scratch_concat.xtc -cat
        INPUT_XTC="scratch_concat.xtc"
    else
        INPUT_XTC=\$(ls inputs/*.xtc | head -n 1)
    fi

    # 2. Generate index directly from system topology
    echo "q" | gmx make_ndx -f ${tpr} -o index.ndx

    # 3. Build strict [ COMPLEX ] group (Protein + Ligand + Cofactors, stripped of water and ions)
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action ensure_complex

    # 4. PBC Step 1: Make all molecules whole across periodic boundaries
    echo "0" | gmx trjconv -s ${tpr} -f \${INPUT_XTC} -o pbc_whole.xtc -pbc whole

    # 5. PBC Step 2: Center on Protein, make molecules compact, and extract full COMPLEX
    echo -e "${pbc_center}\\n${pbc_extract}" | gmx trjconv -s ${tpr} -f pbc_whole.xtc -n index.ndx -o pbc_centered.xtc -center -pbc mol -ur compact

    # 6. PBC Step 3: Rotational and translational fitting on Protein backbone to eliminate drift
    echo -e "Protein\\n${pbc_extract}" | gmx trjconv -s ${tpr} -f pbc_centered.xtc -n index.ndx -o trajectory_clean.xtc -fit rot+trans

    # 7. Extract reference 3D structure of the dry COMPLEX at t=0 as standard GROMACS GRO
    echo "${pbc_extract}" | gmx trjconv -s ${tpr} -f trajectory_clean.xtc -n index.ndx -o complex_ref.gro -dump 0

    # Immediately delete intermediate heavy solvent trajectories from scratch
    rm -f pbc_whole.xtc pbc_centered.xtc scratch_concat.xtc sorted_chunks.list
    """

    stub:
    """
    touch trajectory_clean.xtc
    touch complex_ref.gro
    """
}

