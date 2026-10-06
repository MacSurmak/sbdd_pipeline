/*
 * workflows/screen.nf
 * Decoupled parallel docking workflow.
 * Consumes pre-prepared receptors and pre-prepared ligands directly from DB plan.
 */

include { DOCK_LIBRARY      } from "../subworkflows/dock_library"
include { SUMMARIZE_DOCKING } from "../modules/reporting/summarize"
include { SPLIT_SDF         } from "../modules/common/helpers"

workflow SCREEN_WF {
    if (!params.dock_plan) { error "FATAL: Parameter --dock_plan (JSON) is required." }

    def planFile = file(params.dock_plan)
    def jsonSlurper = new groovy.json.JsonSlurper()
    def planData = jsonSlurper.parseText(planFile.text)

    // 1. Receptors
    ch_targets = Channel.fromList(planData.receptors).map { rec ->
        def target_id = rec.name
        def rec_mae = file(rec.path_mae)
        def rec_pdb = file(rec.path_pdb)

        if (!rec_mae.exists()) error "FATAL: Receptor MAE file not found: ${rec.path_mae}"
        if (!rec_pdb.exists()) error "FATAL: Receptor PDB file not found: ${rec.path_pdb}"

        def ref_sdf = file(rec.ref_sdf)
        def ref_mae = file(rec.ref_mae)

        if (!ref_sdf.exists()) error "FATAL: Reference ligand SDF not found for ${target_id}: ${rec.ref_sdf}"
        if (!ref_mae.exists()) error "FATAL: Reference ligand MAE not found for ${target_id}: ${rec.ref_mae}"

        def meta = [
            id: target_id,
            frame: 1,
            project: params.project_name,
            native_ligand: rec.native_ligand ?: ""
        ]

        [ meta, rec_pdb.toAbsolutePath(), rec_mae.toAbsolutePath(), ref_sdf.toAbsolutePath(), ref_mae.toAbsolutePath() ]
    }

    ch_rec_pdb = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, pdb] }
    ch_rec_mae = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, mae] }
    ch_ref_sdf = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, rsdf] }
    ch_ref_mae = ch_targets.map { meta, pdb, mae, rsdf, rmae -> [meta, rmae] }

    // 2. Pre-prepared Ligands
    ch_ready_ligands = Channel.fromList(planData.ligands).map { lig ->
        def sdf_f = file(lig.path_sdf)
        def mae_f = file(lig.path_mae)
        if (!sdf_f.exists()) error "FATAL: Prepared ligand SDF not found: ${lig.path_sdf}"
        if (!mae_f.exists()) error "FATAL: Prepared ligand MAE not found: ${lig.path_mae}"

        def meta = [
            id: lig.name,
            project: params.project_name,
            receptor: lig.receptor ?: "",
            ligand_title: lig.ligand_title ?: lig.name
        ]
        [ meta, sdf_f.toAbsolutePath(), mae_f.toAbsolutePath() ]
    }

    // Chunk prepared SDF for cluster GPU parallelization
    ch_sdf_for_chunk = ch_ready_ligands.map { meta, sdf, mae -> [ meta, sdf ] }
    SPLIT_SDF(ch_sdf_for_chunk, params.batch_size ?: 10)

    ch_chunks_sdf = SPLIT_SDF.out.chunks.flatMap { meta, files ->
        def fileList = files instanceof List ? files : (files ? [files] : [])
        fileList.withIndex().collect { f, idx ->
            def chunk_meta = meta + [chunk_idx: idx]
            [ chunk_meta, f, idx ]
        }
    }

    ch_chunks_mae = ch_ready_ligands.map { meta, sdf, mae ->
        def chunk_meta = meta + [chunk_idx: 0]
        [ chunk_meta, mae ]
    }

    // 3. Parallel Docking
    DOCK_LIBRARY(
        ch_rec_pdb,
        ch_ref_sdf,
        ch_chunks_sdf,
        ch_rec_mae,
        ch_ref_mae,
        ch_chunks_mae,
        params.dock_engine,
        params.exhaustiveness,
        params.cnn_scoring,
        params.num_modes,
        params.autobox_add,
        (planData.mode == "self")
    )

    // 4. Summarize results
    def bench_ligands_str = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""
    SUMMARIZE_DOCKING(
        params.project_name,
        DOCK_LIBRARY.out.summary_files,
        params.native_ligand ?: "",
        bench_ligands_str
    )
}
