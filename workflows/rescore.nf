/*
 * workflows/rescore.nf
 * Standalone Prime MM-GBSA and Consensus Rescoring Workflow
 */

include { BUILD_FUNNEL_WHITELIST }  from '../modules/reporting/build_funnel'
include { FILTER_PV }               from '../modules/schrodinger/filter_pv'
include { PRIME_MMGBSA }            from '../modules/schrodinger/prime_mmgbsa'
include { SUMMARIZE_DOCKING }       from '../modules/reporting/summarize'
include { MERGE_CSV }               from '../modules/common/helpers'

workflow RESCORE_WF {
    if (!params.rescore_plan) {
        log.error "FATAL: Parameter --rescore_plan (JSON) is required."
        return
    }

    def planFile = file(params.rescore_plan)
    if (!planFile.exists()) {
        log.error "FATAL: Rescore plan file not found: ${params.rescore_plan}"
        return
    }

    def jsonSlurper = new groovy.json.JsonSlurper()
    def planData = jsonSlurper.parseText(planFile.text)

    ch_tasks = Channel.fromList(planData.tasks ?: []).map { task ->
        def pv_file = file(task.pv_path)
        if (!pv_file.exists()) error "FATAL: Pose Viewer file not found: ${task.pv_path}"

        def meta = [
            id: task.name,
            frame: task.frame ?: 1,
            project: params.project_name
        ]
        [ meta, pv_file.toAbsolutePath() ]
    }

    def prime_job = params.prime_job_type ?: "ENERGY"
    def mem_enabled = params.membrane ?: false

    if (planData.use_whitelist && planData.whitelist_path) {
        ch_wl = Channel.fromPath(planData.whitelist_path)
        ch_filter_in = ch_tasks.combine(ch_wl).map { meta, pv, wl -> [ meta, pv ] }

        FILTER_PV(ch_filter_in, ch_wl, params.mmgbsa_batch_size ?: 10)

        ch_mmgbsa_in = FILTER_PV.out.filtered_pv_chunks.flatMap { meta, files ->
            def fileList = files instanceof List ? files : [files]
            fileList.sort { it.name }.withIndex().collect { f, idx ->
                def chunk_meta = [
                    id: "${meta.id}_mmc${idx}",
                    frame_id: meta.id,
                    frame: meta.frame,
                    chunk_idx: idx,
                    project: meta.project
                ]
                [ chunk_meta, f ]
            }
        }
    } else {
        ch_mmgbsa_in = ch_tasks
    }

    // Single unified invocation of PRIME_MMGBSA
    PRIME_MMGBSA(ch_mmgbsa_in, prime_job, mem_enabled)

    if (planData.use_whitelist && planData.whitelist_path) {
        ch_to_merge = PRIME_MMGBSA.out.scores_csv.map { chunk_meta, csv ->
            def frame_meta = [id: chunk_meta.frame_id, frame: chunk_meta.frame, project: chunk_meta.project]
            [ frame_meta, csv ]
        }.groupTuple(by: 0, sort: true)

        MERGE_CSV(ch_to_merge)
        ch_mmgbsa_csvs = MERGE_CSV.out.merged_csv
        ch_mmgbsa_maes = PRIME_MMGBSA.out.out_maegz.map { chunk_meta, mae ->
            def frame_meta = [id: chunk_meta.frame_id, frame: chunk_meta.frame, project: chunk_meta.project]
            [ frame_meta, mae ]
        }.groupTuple(by: 0, sort: true).map { frame_meta, maes -> [ frame_meta, maes[0] ] }
    } else {
        ch_mmgbsa_csvs = PRIME_MMGBSA.out.scores_csv
        ch_mmgbsa_maes = PRIME_MMGBSA.out.out_maegz
    }

    // Generate execution receipt for atomic database synchronization targeting published deliverables
    ch_receipt_entries = ch_mmgbsa_csvs.join(ch_mmgbsa_maes, by: 0).map { meta, csv, mae ->
        def frame_id = meta.frame_id ?: meta.id
        def target_dir = meta.ligand_title ? "${params.outdir}/docking/${frame_id}/${meta.ligand_title}" : "${params.outdir}/docking/${frame_id}"
        "${meta.id}\t${target_dir}/${csv.name}\t${target_dir}/${mae.name}"
    }
    ch_receipt_entries.collectFile(name: 'rescore_receipt.tsv', storeDir: params.outdir, newLine: true)

    // Re-summarize master ranking with deterministic staging from plan
    def bench_str = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""

    if (!planData.dock_summary_files) {
        error "FATAL: 'dock_summary_files' list is missing in rescore_plan.json. Manifest must be strictly deterministic."
    }

    ch_dock_files = Channel.fromList(planData.dock_summary_files).map { path_str ->
        def f = file(path_str)
        if (!f.exists()) {
            error "FATAL: Required docking artifact referenced in plan does not exist: ${path_str}"
        }
        f
    }

    ch_rescore_plan = Channel.fromPath(params.rescore_plan, checkIfExists: true)

    ch_all_summary = ch_mmgbsa_csvs.map { meta, csv -> csv }
        .mix(ch_dock_files)
        .collect()

    SUMMARIZE_DOCKING(
        params.project_name,
        ch_rescore_plan,
        ch_all_summary,
        params.native_ligand ?: "",
        bench_str
    )
}
