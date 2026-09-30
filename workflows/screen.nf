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
        
        // Infer the generated files based on the primary path
        def basePath = file(rec.path).parent
        
        // It could be a MAE or PDB. We know the preparation step generates both.
        def baseName = rec.name
        // Because of auto-suffix _ready from DB, we strip it to find original prepared files, 
        // OR we just use standard naming convention. 
        // In Step 2 we forced DB path to be the .mae or .pdb file.
        def rec_mae = file(rec.path)
        if (!rec_mae.name.endsWith(".mae")) {
            rec_mae = file(rec.path.replace(".pdb", ".mae"))
        }
        def rec_pdb = file(rec.path.replace(".mae", ".pdb").replace(".maegz", ".pdb"))
        
        // Re-construct native ligand path. Since name format is {original_name}_ready,
        def orig_name = rec.name.replaceAll("_ready\$", "")
        def ref_sdf = file("${basePath}/${orig_name}_ref_native_prepared.sdf")
        if (!ref_sdf.exists()) ref_sdf = file("${basePath}/${orig_name}_ref_native.sdf")
        def ref_mae = file("${basePath}/${orig_name}_ref_native_prepared.mae")
        if (!ref_mae.exists()) ref_mae = file("${basePath}/${orig_name}_ref_native.mae")

        if (!ref_sdf.exists() || !ref_mae.exists()) {
            // Provide dummy if absolutely no native ligand
            ref_sdf = file("NO_FILE")
            ref_mae = file("NO_FILE")
        }

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
        [ [id: lig.name, project: params.project_name], file(lig.path).toAbsolutePath() ] 
    }

    PREPARE_LIGANDS(
        ch_lib,
        ch_ref_mae,
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
