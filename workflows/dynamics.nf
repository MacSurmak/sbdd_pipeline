/*
 * workflows/dynamics.nf
 * Standalone MD system setup, chunked simulation, and trajectory analysis
 */

include { RUN_CHUNKED_MD   } from '../subworkflows/run_chunked_md'

workflow DYNAMICS_WF {
    if (!params.md_protocol) { error "FATAL: Parameter --md_protocol is required. Did you run 'sbdd md-init'?" }
    if (!params.md_plan) { error "FATAL: Parameter --md_plan must be provided for dynamics." }

    def planFile = file(params.md_plan)
    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    ch_systems_to_simulate = Channel.fromList(planData).map { target ->
        def meta = [
            id: target.id,
            project: params.project_name,
            replicates: target.reps_to_run
        ]
        def gro = file(target.gro).toAbsolutePath()
        def top = file(target.top).toAbsolutePath()
        def ndx = file(target.ndx).toAbsolutePath()
        def posre_files = file("${target.top_dir}/posre_*.itp").collect { it.toAbsolutePath() }

        [ meta, gro, top, ndx, posre_files ]
    }

    ch_md_protocol = Channel.fromPath("${params.md_protocol}/*", checkIfExists: true).collect()

    // Equilibration & Production using Protocol Playlist
    RUN_CHUNKED_MD(
        ch_systems_to_simulate,
        ch_md_protocol,
        params.md_total_ns,
        params.md_chunk_ns,
        params.md_replicates as int
    )
}
