/*
 * workflows/screen.nf
 * Standalone decoupled distributed docking workflow.
 * Reads directly from manifest targets and prepared/ directory.
 */

include { PREPARE_LIGANDS   } from "../subworkflows/prepare_ligands"
include { DOCK_LIBRARY      } from "../subworkflows/dock_library"
include { SUMMARIZE_DOCKING } from "../modules/reporting/summarize"

workflow SCREEN_WF {
    if (!params.manifest) { error "FATAL: Parameter --manifest is required." }
    if (!params.library) { error "FATAL: Parameter --library is required for screening." }

    def manifestFile = file(params.manifest)
    def yaml = new org.yaml.snakeyaml.Yaml()
    def manifestData = yaml.load(manifestFile.text)

    // 1. Load pre-prepared structures from the prepared/ folder
    ch_targets = Channel.fromList(manifestData.targets).map { target ->
        def target_id = target.id
        
        def rec_pdb = file("${params.outdir}/prepared/${target_id}_prepared.pdb")
        def rec_mae = file("${params.outdir}/prepared/${target_id}_prepared.mae")
        def ref_sdf = file("${params.outdir}/prepared/${target_id}_ref_native.sdf")
        def ref_mae = file("${params.outdir}/prepared/${target_id}_ref_native.mae")
        
        if (!rec_pdb.exists() || !rec_mae.exists()) {
            error "FATAL: Prepared structures not found for target ${target_id} in prepared/ directory. Did you run 'sbdd prep'?"
        }

        def meta = [
            id: target_id,
            frame: 1, // Currently assuming single frame for decoupled screen
            project: params.project_name
        ]
        
        [ meta, rec_pdb.toAbsolutePath(), rec_mae.toAbsolutePath(), ref_sdf.toAbsolutePath(), ref_mae.toAbsolutePath() ]
    }

    // Unpack channel tuples for explicit arguments
    ch_rec_pdb = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, pdb] }
    ch_rec_mae = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, mae] }
    ch_ref_sdf = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, rsdf] }
    ch_ref_mae = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, rmae] }

    // 2. Prepare Library
    ch_lib = Channel.fromPath(file(params.library).toAbsolutePath(), checkIfExists: true)
        .map { lib -> [ [id: params.project_name, project: params.project_name], lib ] }

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
