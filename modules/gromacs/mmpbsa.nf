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
    path(templates_dir)
    val(start_ns)
    val(sample_interval_ps)
    val(temp)
    val(salt_conc)
    val(is_membrane)
    val(decomp)
    val(method)
    val(gb_model)
    val(emem)
    val(indi)
    val(fillratio)

    output:
    tuple val(meta), path("thermodynamics_mmpbsa_pb.csv"), path("thermodynamics_mmpbsa_pb.dat"), emit: pb_results
    tuple val(meta), path("decomposition_mmpbsa_gb.csv"), path("thermodynamics_mmpbsa_gb.dat"), emit: gb_decomp_results, optional: true

    script:
    if (is_membrane.toString() == 'true' && !meta.membrane_thickness) {
        error "FATAL: Membrane is enabled for ${meta.id}, but meta.membrane_thickness is missing!"
    }
    def mem_arg = (is_membrane.toString() == 'true') ? "--membrane --membrane_thickness ${meta.membrane_thickness}" : ""
    def dec_arg = (decomp.toString() == 'true') ? "--decomp" : ""
    """
    python3 ${projectDir}/bin/wrap_gromacs.py mmpbsa \\
        --clean_xtc ${clean_xtc} \\
        --ref_pdb ${ref_pdb} \\
        --topol ${topol_top} \\
        --ndx ${index_ndx} \\
        --ligand "${meta.ligand}" \\
        --templates_dir "${templates_dir}" \\
        --start_ns ${start_ns} \\
        --sample_interval_ps ${sample_interval_ps} \\
        --temp ${temp} \\
        --salt_conc ${salt_conc} \\
        --method "${method}" \\
        --gb_model ${gb_model} \\
        --emem ${emem} \\
        --indi ${indi} \\
        --fillratio ${fillratio} \\
        --cpus ${task.cpus} \\
        ${mem_arg} \\
        ${dec_arg}
    """
}

