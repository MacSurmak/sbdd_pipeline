/*
 * modules/gromacs/mmpbsa.nf
 * Free Energy Calculation (Membrane PB) and Per-Residue Decomposition (GB) via gmx_MMPBSA
 * Accurately reproducing original physical parameters and thermodynamic models.
 */

process GMX_MMPBSA {
    tag "${meta.id}_rep${meta.rep}"
    label 'process_high'

    input:
    tuple val(meta), path(clean_xtc), path(ref_pdb), path(topol_top), path(index_ndx)
    val(start_ns)
    val(step_frames)
    val(temp)
    val(salt_conc)
    val(is_membrane)
    val(membrane_thickness)
    val(decomp)

    output:
    tuple val(meta), path("thermodynamics_mmpbsa_pb.csv"), path("thermodynamics_mmpbsa_pb.dat"), emit: pb_results
    tuple val(meta), path("decomposition_mmpbsa_gb.csv"), path("thermodynamics_mmpbsa_gb.dat"), emit: gb_decomp_results, optional: true

    script:
    def mem_arg = (is_membrane.toString() == 'true') ? "--membrane --membrane_thickness ${membrane_thickness}" : ""
    def dec_arg = (decomp.toString() == 'true') ? "--decomp" : ""
    """
    python3 ${projectDir}/bin/wrap_gromacs.py mmpbsa \\
        --clean_xtc ${clean_xtc} \\
        --ref_pdb ${ref_pdb} \\
        --topol ${topol_top} \\
        --ndx ${index_ndx} \\
        --ligand "${meta.ligand}" \\
        --start_ns ${start_ns} \\
        --step_frames ${step_frames} \\
        --temp ${temp} \\
        --salt_conc ${salt_conc} \\
        --cpus ${task.cpus} \\
        ${mem_arg} \\
        ${dec_arg}
    """

    stub:
    """
    touch thermodynamics_mmpbsa_pb.csv
    touch thermodynamics_mmpbsa_pb.dat
    touch decomposition_mmpbsa_gb.csv
    touch thermodynamics_mmpbsa_gb.dat
    """
}

