/*
 * subworkflows/analyze_md_traj.nf
 * Analyzes clean trajectories (RMSD, RMSF, Contacts) and runs gmx_MMPBSA
 */

include { MD_METRICS } from '../modules/gromacs/md_metrics'
include { GMX_MMPBSA } from '../modules/gromacs/mmpbsa'

workflow ANALYZE_MD_TRAJ {
    take:
    ch_traj_clean      // channel: [ val(meta), path(clean_xtc), path(ref_pdb), path(topol), path(ndx) ]
    start_ns           // val: 20.0
    step_frames        // val: 40
    temperature        // val: 310.15
    salt_conc          // val: 0.15
    is_membrane        // val: true/false
    membrane_thickness // val: 28.0
    per_res_decomp     // val: true

    main:
    // 1. Calculate structural metrics
    ch_metrics_in = ch_traj_clean.map { meta, clean_xtc, ref_pdb, topol, ndx ->
        [ meta, clean_xtc, ref_pdb ]
    }
    MD_METRICS(ch_metrics_in, params.metric_heavy, params.metric_hyd, params.metric_pol)

    // 2. Run MM-PBSA free energy calculation
    GMX_MMPBSA(
        ch_traj_clean,
        start_ns,
        step_frames,
        temperature,
        salt_conc,
        is_membrane,
        membrane_thickness,
        per_res_decomp
    )

    emit:
    metrics_tables = MD_METRICS.out.tables
    metrics_json   = MD_METRICS.out.summary_json
    pb_results     = GMX_MMPBSA.out.pb_results
    gb_decomp      = GMX_MMPBSA.out.gb_decomp_results
}

