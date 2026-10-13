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
    val(water_names)
    val(ion_names)
    val(lipid_names)

    output:
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.pdb"), emit: clean_traj

    script:
    """
    python3 ${projectDir}/bin/wrap_gromacs.py trjconv \\
        --xtcs inputs/*.xtc \\
        --tpr ${tpr} \\
        --pbc_center ${pbc_center} \\
        --pbc_extract ${pbc_extract} \\
        --out_xtc trajectory_clean.xtc \\
        --out_pdb complex_ref.pdb \\
        --water_names "${water_names}" \\
        --ion_names "${ion_names}" \\
        --lipid_names "${lipid_names}"
    """

    stub:
    """
    touch trajectory_clean.xtc
    touch complex_ref.pdb
    """
}

