/*
 * modules/gromacs/trjconv.nf
 * PBC correction, rotational/translational fitting, and water/ion stripping
 */

process TRJCONV {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(xtc_inputs, stageAs: 'inputs/*'), path(tpr), path(dry_top)
    val(pbc_center)
    val(pbc_extract)
    val(water_names)
    val(ion_names)
    val(lipid_names)

    output:
    tuple val(meta), path("trajectory_clean.xtc"), path("complex_ref.pdb"), path("complex_clean.ndx"), emit: clean_traj

    script:
    def lig_arg = meta.ligand ? "--ligand ${meta.ligand}" : ""
    def dry_arg = (dry_top && dry_top.name != 'NO_DRY_TOP') ? "--dry_top ${dry_top}" : ""
    """
    python3 ${projectDir}/bin/wrap_gromacs.py trjconv \\
        --xtcs inputs/*.xtc \\
        --tpr ${tpr} \\
        --pbc_center ${pbc_center} \\
        --pbc_extract ${pbc_extract} \\
        --out_xtc trajectory_clean.xtc \\
        --out_pdb complex_ref.pdb \\
        --out_ndx complex_clean.ndx \\
        ${lig_arg} \\
        ${dry_arg} \\
        --water_names "${water_names}" \\
        --ion_names "${ion_names}" \\
        --lipid_names "${lipid_names}"
    """
}

