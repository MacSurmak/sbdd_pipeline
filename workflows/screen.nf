/*
 * workflows/screen.nf
 * Standalone decoupled distributed docking workflow.
 * Reads directly from manifest targets and prepared/ directory.
 */

include { PREPARE_LIGANDS   } from "../subworkflows/prepare_ligands"
include { DOCK_LIBRARY      } from "../subworkflows/dock_library"
include { SUMMARIZE_DOCKING } from "../modules/reporting/summarize"

workflow SCREEN_WF {
    if (!params.dock_plan) { error "FATAL: Parameter --dock_plan (JSON) is required." }

    def planFile = file(params.dock_plan)
    def jsonSlurper = new groovy.json.JsonSlurper()
    def planData = jsonSlurper.parseText(planFile.text)

    // 1. Parse Receptors (They are already strictly verified as 'ready' with absolute paths by CLI)
    ch_targets = Channel.fromList(planData.receptors).map { rec ->
        def target_id = rec.name
        def rec_path = file(rec.path)
        def rec_mae = rec_path
        if (!rec_mae.name.endsWith(".mae")) {
            rec_mae = file(rec.path.replace(".pdb", ".mae"))
        }
        def rec_pdb = file(rec.path.replace(".mae", ".pdb").replace(".maegz", ".pdb"))

        def ref_path = (rec.ref_pose_path && rec.ref_pose_path != "NO_FILE") ? file(rec.ref_pose_path) : file("NO_FILE")
        def ref_sdf = ref_path
        def ref_mae = ref_path

        def meta = [
            id: target_id,
            frame: 1, 
            project: params.project_name
        ]
        
        [ meta, rec_pdb.toAbsolutePath(), rec_mae.toAbsolutePath(), ref_sdf.toAbsolutePath(), ref_mae.toAbsolutePath() ]
    }

    ch_rec_pdb = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, pdb] }
    ch_rec_mae = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, mae] }
    ch_ref_sdf = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, rsdf] }
    ch_ref_mae = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, rmae] }

    // 2. Parse Ligands (Support for multiple libraries simultaneously)
    ch_lib = Channel.fromList(planData.ligands).map { lig -> 
        [ [id: lig.name, project: params.project_name, status: lig.status, is_native: lig.is_native], file(lig.path).toAbsolutePath() ] 
    }

    PREPARE_LIGANDS(
        ch_lib,
        params.dock_engine,
        params.batch_size,
        params.ligand_ph,
        params.ligand_ph_tol,
        params.ligprep_ff,
        params.max_states
    )

    // 3. Parallel Docking
    DOCK_LIBRARY(
        ch_rec_pdb,
        ch_ref_sdf,
        PREPARE_LIGANDS.out.chunks,
        ch_rec_mae,
        ch_ref_mae,
        PREPARE_LIGANDS.out.prepared_ligands,
        params.dock_engine,
        params.exhaustiveness,
        params.cnn_scoring,
        params.num_modes,
        params.autobox_add,
        params.top_k_frames,
        params.percentile_cutoff,
        params.membrane
    )

    // 4. Summarize results
    def bench_ligands_str = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""
    SUMMARIZE_DOCKING(
        params.project_name,
        DOCK_LIBRARY.out.summary_files,
        "", // Default native ligand extraction handled gracefully
        bench_ligands_str
    )
}
