/*
 * workflows/dynamics.nf
 * Standalone MD system setup, chunked simulation, and trajectory analysis
 */

include { STAGE_INPUTS    } from '../modules/common/stage_inputs'
include { SETUP_MD_SYSTEM } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD  } from '../subworkflows/run_chunked_md'
include { ANALYZE_MD_TRAJ } from '../subworkflows/analyze_md_traj'

workflow DYNAMICS_WF {
    def mdp_minim = file("${projectDir}/assets/mdp/minim.mdp")
    def mdp_nvt   = file("${projectDir}/assets/mdp/nvt.mdp")
    def mdp_npt   = file("${projectDir}/assets/mdp/npt.mdp")
    def mdp_md    = file("${projectDir}/assets/mdp/md.mdp")

    // Automatic resolution of inputs
    def prepared_pdb = file("${params.outdir}/receptor/${params.project_name}_prepared.pdb")
    def docked_sdf   = file("${params.outdir}/docking/${params.project_name}_docked_all.sdf")
    def ref_sdf      = file("${params.outdir}/receptor/${params.project_name}_ref_native.sdf")

    def rec_file     = prepared_pdb.exists() ? prepared_pdb : file(params.receptor)
    def dock_file    = docked_sdf.exists() ? docked_sdf : (params.docked_poses ? file(params.docked_poses) : ref_sdf)

    STAGE_INPUTS(rec_file, dock_file, ref_sdf)

    // Parse target hit list (e.g. "native,PPZ_1H_Ph_BnMe" or single hit)
    def hits = params.md_hits ? params.md_hits.split(",").collect { it.trim() } : ["native"]

    ch_targets = Channel.fromList(hits)
        .combine(STAGE_INPUTS.out.staged_library)
        .combine(STAGE_INPUTS.out.staged_receptor)
        .map { hit_name, staged_dock, staged_rec ->
            def target_id = hit_name.replaceAll("[^a-zA-Z0-9_-]", "_")
            [ [id: target_id], staged_dock, hit_name, staged_rec ]
        }

    // 1. Build systems
    SETUP_MD_SYSTEM(
        ch_targets,
        params.membrane ?: false,
        params.upper_lipids ?: "POPC:POPE:CHL1",
        params.upper_ratio ?: "5:4:1",
        params.lower_lipids ?: "POPC:POPE:POPS:CHL1",
        params.lower_ratio ?: "4:4:1:1",
        10.0,
        15.0,
        params.md_salt_conc ?: 0.15,
        "opc"
    )

    // 2. Equilibration & chunked production MD on Slurm GPU
    RUN_CHUNKED_MD(
        SETUP_MD_SYSTEM.out.system_ready,
        mdp_minim,
        mdp_nvt,
        mdp_npt,
        mdp_md,
        params.membrane ?: false,
        params.md_temperature ?: 310.15,
        params.md_total_ns ?: 60,
        params.md_chunk_ns ?: 10,
        params.md_replicates ?: 3
    )

    // 3. Trajectory analysis & MM-PBSA
    ANALYZE_MD_TRAJ(
        RUN_CHUNKED_MD.out.clean_trajectory,
        params.mmpbsa_cores ?: 16,
        params.mmpbsa_start_ns ?: 20.0,
        params.mmpbsa_step_frames ?: 40,
        params.md_temperature ?: 310.15,
        true
    )
}
