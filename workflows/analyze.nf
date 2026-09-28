/*
 * workflows/analyze.nf
 * Independent trajectory analysis workflow (PBC fixing, RMSD, MMPBSA)
 */

include { TRJCONV         } from '../modules/gromacs/trjconv'
include { ANALYZE_MD_TRAJ } from '../subworkflows/analyze_md_traj'

workflow ANALYZE_WF {
    if (!params.traj || !params.tpr || !params.topol || !params.ndx) {
        error "Analysis workflow requires explicitly provided --traj, --tpr, --topol, and --ndx parameters."
    }

    ch_traj  = Channel.fromPath(params.traj, checkIfExists: true)
    ch_tpr   = Channel.fromPath(params.tpr, checkIfExists: true)
    ch_topol = Channel.fromPath(params.topol, checkIfExists: true)
    ch_ndx   = Channel.fromPath(params.ndx, checkIfExists: true)

    // Dummy meta block since we operate on single trajectories in standalone mode
    def meta = [id: params.project_name ?: "md_analysis", rep: 1]

    ch_trjconv_in = Channel.value(meta)
        .combine(ch_traj)
        .combine(ch_tpr)
        .combine(ch_topol)
        .combine(ch_ndx)
        .map { m, traj, tpr, top, ndx -> [m, traj, tpr, top, ndx] }

    TRJCONV(
        ch_trjconv_in,
        params.pbc_center,
        params.pbc_extract
    )

    ANALYZE_MD_TRAJ(
        TRJCONV.out.clean_traj,
        params.mmpbsa_start_ns,
        params.mmpbsa_step_frames,
        params.md_temperature,
        params.md_salt_conc,
        params.membrane,
        params.membrane_thickness,
        params.per_residue_decomp
    )
}
