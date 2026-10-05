/*
 * workflows/rescore.nf
 * Standalone Prime MM-GBSA and Consensus Rescoring Workflow
 */

include { BUILD_FUNNEL_WHITELIST }  from '../modules/reporting/build_funnel'
include { FILTER_PV }               from '../modules/schrodinger/filter_pv'
include { PRIME_MMGBSA }            from '../modules/schrodinger/prime_mmgbsa'
include { SUMMARIZE_DOCKING }       from '../modules/reporting/summarize'

process MERGE_MMGBSA_CHUNKS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(csv_files, stageAs: 'mmgbsa_?/*')

    output:
    tuple val(meta), path("${meta.id}_mmgbsa-out.csv"), emit: scores_csv

    script:
    """
    python3 ${projectDir}/bin/merge_csv.py --inputs ${csv_files} --output ${meta.id}_mmgbsa-out.csv
    """

    stub:
    """
    touch ${meta.id}_mmgbsa-out.csv
    """
}

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

        MERGE_MMGBSA_CHUNKS(ch_to_merge)
        ch_all_csvs = MERGE_MMGBSA_CHUNKS.out.scores_csv.map { meta, csv -> csv }.collect()
    } else {
        ch_all_csvs = PRIME_MMGBSA.out.scores_csv.map { meta, csv -> csv }.collect()
    }

    // Re-summarize master ranking including Prime MM-GBSA
    def bench_str = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""

    SUMMARIZE_DOCKING(
        params.project_name,
        ch_all_csvs,
        params.native_ligand ?: "",
        bench_str
    )
}
