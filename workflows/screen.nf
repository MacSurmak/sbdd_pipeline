/*
 * workflows/screen.nf
 * Decoupled parallel docking workflow.
 * Consumes pre-prepared receptors and pre-prepared ligands directly from DB plan.
 */

include { DOCK_LIBRARY       } from "../subworkflows/dock_library"
include { SUMMARIZE_DOCKING  } from "../modules/reporting/summarize"
include { BATCH_DOCK_LIGANDS } from "../modules/common/helpers"

workflow SCREEN_WF {
    if (!params.dock_plan) { error "FATAL: Parameter --dock_plan (JSON) is required." }
    if (params.glide_batch_size == null || params.gnina_batch_size == null) {
        error "FATAL: Parameters --glide_batch_size and --gnina_batch_size are required for docking."
    }

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

    // 2. Batch Pre-prepared Ligands into independent chunks for Glide (.maegz) and GNINA (.sdf)
    ch_dock_plan_file = Channel.fromPath(params.dock_plan, checkIfExists: true)
    BATCH_DOCK_LIGANDS(
        params.project_name,
        ch_dock_plan_file,
        params.glide_batch_size,
        params.gnina_batch_size
    )

    ch_glide_manifest = BATCH_DOCK_LIGANDS.out.manifest.flatMap { mf ->
        def parsed = new groovy.json.JsonSlurperClassic().parseText(mf.text)
        parsed.glide_chunks ?: []
    }
    ch_gnina_manifest = BATCH_DOCK_LIGANDS.out.manifest.flatMap { mf ->
        def parsed = new groovy.json.JsonSlurperClassic().parseText(mf.text)
        parsed.gnina_chunks ?: []
    }

    ch_sdf_map = BATCH_DOCK_LIGANDS.out.sdf_chunks.flatten().map { f -> [ f.name, f ] }
    ch_mae_map = BATCH_DOCK_LIGANDS.out.mae_chunks.flatten().map { f -> [ f.name, f ] }

    ch_chunks_sdf = ch_gnina_manifest
        .map { item -> [ item.file, item ] }
        .join(ch_sdf_map, by: 0)
        .map { sdf_name, item, sdf_file ->
            def chunk_meta = [
                id: "c${item.chunk_idx}",
                chunk_idx: item.chunk_idx,
                receptor: item.receptor ?: "",
                project: params.project_name
            ]
            [ chunk_meta, sdf_file, item.chunk_idx ]
        }

    ch_chunks_mae = ch_glide_manifest
        .map { item -> [ item.file, item ] }
        .join(ch_mae_map, by: 0)
        .map { mae_name, item, mae_file ->
            def chunk_meta = [
                id: "c${item.chunk_idx}",
                chunk_idx: item.chunk_idx,
                receptor: item.receptor ?: "",
                project: params.project_name
            ]
            [ chunk_meta, mae_file ]
        }

    // Validate Fail-Fast parameters before dispatching tasks
    if (params.dock_engine in ["glide", "both"]) {
        if (!params.inner_box || !params.outer_box || !params.glide_forcefield || !params.glide_precision || !params.poses_per_lig || params.glide_postdock == null || !params.glide_recep_vscale || !params.glide_lig_vscale || !params.glide_recep_ccut || !params.glide_lig_ccut) {
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
        params.glide_postdock,
        params.glide_recep_vscale,
        params.glide_lig_vscale,
        params.glide_recep_ccut,
        params.glide_lig_ccut,
        (planData.mode == "self"),
        params.outdir,
        params.dock_run_name
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
