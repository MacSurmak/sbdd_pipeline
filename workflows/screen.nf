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
    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

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

        def ion_sdfs = rec.ion_sdfs ? rec.ion_sdfs.collect { file(it) } : []
        def ion_maes = rec.ion_maes ? rec.ion_maes.collect { file(it) } : []
        def cof_sdfs = rec.cofactor_sdfs ? rec.cofactor_sdfs.collect { file(it) } : []
        def cof_maes = rec.cofactor_maes ? rec.cofactor_maes.collect { file(it) } : []

        def meta = [
            id: target_id,
            frame: 1,
            project: params.project_name,
            native_ligand: rec.native_ligand ?: "",
            ion_resnames: rec.ion_resnames ?: [],
            cofactor_resnames: rec.cofactor_resnames ?: []
        ]

        [ meta, rec_pdb.toAbsolutePath(), rec_mae.toAbsolutePath(), ref_sdf.toAbsolutePath(), ref_mae.toAbsolutePath(), ion_sdfs*.toAbsolutePath(), ion_maes*.toAbsolutePath(), cof_sdfs*.toAbsolutePath(), cof_maes*.toAbsolutePath() ]
    }

    ch_rec_pdb  = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, pdb] }
    ch_rec_mae  = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, mae] }
    ch_ref_sdf  = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, rsdf] }
    ch_ref_mae  = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, rmae] }
    ch_ion_sdfs = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, isdf] }
    ch_ion_maes = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, imae] }
    ch_cof_sdfs = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, csdf] }
    ch_cof_maes = ch_targets.map { meta, pdb, mae, rsdf, rmae, isdf, imae, csdf, cmae -> [meta, cmae] }

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

    // Validate Fail-Fast parameters before dispatching tasks
    if (params.dock_engine in ["glide", "both"]) {
        if (!params.inner_box || !params.outer_box || !params.glide_forcefield || !params.glide_precision || !params.poses_per_lig || !params.glide_recep_vscale || !params.glide_lig_vscale || !params.glide_recep_ccut || !params.glide_lig_ccut) {
            error "FATAL: Glide docking parameters missing in configuration! Check glide.yaml."
        }
    }
    if (params.dock_engine in ["gnina", "both"]) {
        if (!params.exhaustiveness || !params.cnn_scoring || !params.num_modes || !params.autobox_add) {
            error "FATAL: GNINA docking parameters missing in configuration! Check gnina.yaml."
        }
    }

    // 3. Parallel Docking
    DOCK_LIBRARY(
        ch_rec_pdb,
        ch_ref_sdf,
        ch_chunks_sdf,
        ch_rec_mae,
        ch_ref_mae,
        ch_chunks_mae,
        ch_ion_sdfs,
        ch_ion_maes,
        ch_cof_sdfs,
        ch_cof_maes,
        params.dock_engine,
        params.exhaustiveness,
        params.cnn_scoring,
        params.num_modes,
        params.autobox_add,
        params.inner_box,
        params.outer_box,
        params.glide_forcefield,
        params.glide_precision,
        params.poses_per_lig,
        params.cnn_model ?: "",
        params.glide_recep_vscale,
        params.glide_lig_vscale,
        params.glide_recep_ccut,
        params.glide_lig_ccut,
        (planData.mode == "self"),
        params.outdir
    )

    // 4. Generate structured execution receipt
    DOCK_LIBRARY.out.receipt_entries
        .collectFile(name: 'dock_receipt.tsv', storeDir: params.outdir, newLine: true)

    // 5. Summarize results
    def bench_ligands_str = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""
    ch_dock_plan = Channel.fromPath(params.dock_plan, checkIfExists: true)

    SUMMARIZE_DOCKING(
        params.project_name,
        ch_dock_plan,
        DOCK_LIBRARY.out.summary_files,
        params.native_ligand ?: "",
        bench_ligands_str
    )
}
