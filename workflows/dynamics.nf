/*
 * workflows/dynamics.nf
 * Standalone MD system setup, chunked simulation, and trajectory analysis
 */

include { SETUP_MD_SYSTEM } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD  } from '../subworkflows/run_chunked_md'
include { ANALYZE_MD_TRAJ } from '../subworkflows/analyze_md_traj'
include { STAGE_ASSETS    } from '../modules/common/stage_assets'

workflow DYNAMICS_WF {
    // Automatic resolution of inputs (supporting both single model and multi-frame ensemble naming)
    def prepared_pdb_std = file("${params.outdir}/receptor/${params.project_name}_prepared.pdb")
    def prepared_pdb_f1  = file("${params.outdir}/receptor/${params.project_name}_f1_prepared.pdb")
    def prepared_pdb     = prepared_pdb_std.exists() ? prepared_pdb_std : prepared_pdb_f1

    def docked_sdf   = file("${params.outdir}/docking/${params.project_name}_docked_all.sdf")

    def ref_sdf_std  = file("${params.outdir}/receptor/${params.project_name}_ref_native.sdf")
    def ref_sdf_f1   = file("${params.outdir}/receptor/${params.project_name}_f1_ref_native.sdf")
    def ref_sdf      = ref_sdf_std.exists() ? ref_sdf_std : ref_sdf_f1

    def rec_file     = prepared_pdb.exists() ? prepared_pdb : (params.receptor ? file(params.receptor) : null)
    def dock_file    = docked_sdf.exists() ? docked_sdf : (params.docked_poses ? file(params.docked_poses) : (ref_sdf.exists() ? ref_sdf : null))

    if (!rec_file) { error "No receptor found for dynamics. Provide --receptor or run 'sbdd prepare' first." }
    if (!dock_file) { error "No docked poses found for dynamics. Provide --docked_poses or run 'sbdd screen' first." }

    // Parse target hit list (e.g. "native,PPZ_1H_Ph_BnMe" or single hit)
    def hits = params.md_hits ? params.md_hits.split(",").collect { it.trim() } : ["native"]

    ch_targets = Channel.fromList(hits)
        .map { hit_name ->
            def target_id = hit_name.replaceAll("[^a-zA-Z0-9_-]", "_")
            [ [id: target_id], dock_file, hit_name, rec_file ]
        }

    // 1. Build systems with dynamic padding and water model
    SETUP_MD_SYSTEM(
        ch_targets,
        params.membrane ?: false,
        params.upper_lipids ?: "POPC:POPE:CHL1",
        params.upper_ratio ?: "5:4:1",
        params.lower_lipids ?: "POPC:POPE:POPS:CHL1",
        params.lower_ratio ?: "4:4:1:1",
        params.lipid_padding ?: 10.0,
        params.box_padding ?: 15.0,
        params.md_salt_conc ?: 0.15,
        params.water_model ?: "opc"
    )

    // 2. Equilibration & dynamic atomic production MD
    STAGE_ASSETS(file("${projectDir}/assets/mdp", type: 'dir'))

    RUN_CHUNKED_MD(
        SETUP_MD_SYSTEM.out.system_ready,
        params.membrane ?: false,
        params.md_temperature ?: 310.15,
        params.md_total_ns ?: 60,
        params.md_chunk_ns ?: 10,
        params.md_replicates ? (params.md_replicates as int) : 3,
        params.queue ?: "main",
        params.nvt_short_ps ?: 50,
        params.nvt_long_ps ?: 500,
        params.npt_berendsen_ps ?: 1000,
        params.npt_pr_ps ?: 1000,
        STAGE_ASSETS.out.mdp_templates
    )

    // 3. Trajectory analysis & MM-PBSA
    ANALYZE_MD_TRAJ(
        RUN_CHUNKED_MD.out.clean_trajectory,
        params.mmpbsa_cores ?: 16,
        params.mmpbsa_start_ns ?: 20.0,
        params.mmpbsa_step_frames ?: 40,
        params.md_temperature ?: 310.15,
        params.md_salt_conc ?: 0.15,
        params.membrane ?: false,
        params.membrane_thickness ?: 28.0,
        params.per_residue_decomp ?: true
    )
}
