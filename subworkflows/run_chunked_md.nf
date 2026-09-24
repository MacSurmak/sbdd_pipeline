/*
 * subworkflows/run_chunked_md.nf
 * Replicated production MD simulation with fully dynamic atomic chunk scheduling.
 * Supports any total simulation length and chunk interval (e.g. 500 ns / 1 ns chunks)
 * with zero hardcoded limits.
 */

include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
include { RUN_MD_CHUNK      } from '../modules/gromacs/run_chunk'
include { TRJCAT            } from '../modules/gromacs/trjcat'

workflow RUN_CHUNKED_MD {
    take:
    ch_system        // channel: [ val(meta), path(gro), path(topol), path(ndx), path(posre_files) ]
    ch_mdp_templates // channel: [ path(md.mdp), path(minim.mdp), ... ]
    is_membrane      // val: true/false
    temperature      // val: 310.15
    total_ns         // val: 60
    chunk_ns         // val: 10
    replicates       // val: 3
    slurm_queue      // val: "main"
    nvt_short_ps     // val: 50
    nvt_long_ps      // val: 500
    npt_berendsen_ps // val: 1000
    npt_pr_ps        // val: 1000

    main:
    // 1. Equilibration on Slurm GPU with dynamic step counts
    GMX_EQUILIBRATION(
        ch_system,
        ch_mdp_templates,
        is_membrane,
        temperature,
        nvt_short_ps,
        nvt_long_ps,
        npt_berendsen_ps,
        npt_pr_ps
    )

    def total_c = Math.max(1, Math.round((total_ns as float) / (chunk_ns as float)) as int)
    def num_reps = replicates as int

    // 2. Expand into independent replicates for dynamic production MD
    ch_state = GMX_EQUILIBRATION.out.equilibrated_state
        .join(ch_system.map{ m, gro, top, ndx, pos -> [m.id, top, ndx] }, by: 0)
        .flatMap { target_id, meta, eq_gro, eq_cpt, topol, ndx ->
            (1..num_reps).collect { rep_id ->
                def rep_meta = meta + [rep: rep_id]
                // Send eq_cpt in place of prev_tpr as a dummy for the very first chunk
                [ rep_meta, eq_gro, eq_cpt, eq_cpt, topol, ndx ]
            }
        }

    // 3. Groovy inject() pattern for pure static DAG unrolling
    // Generates a linear chain of tasks exactly total_c times without AST parser errors
    def chunks = (1..total_c).toList()
    def acc = chunks.inject([ch_state, Channel.empty(), Channel.empty()]) { accumulator, c_idx ->
        def current_state = accumulator[0]
        def current_xtcs  = accumulator[1]
        def current_sys   = accumulator[2]

        def step_out = RUN_MD_CHUNK(current_state, ch_mdp_templates, is_membrane, temperature, chunk_ns, c_idx)
        
        def next_state = step_out.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx ->
            [meta, out_gro, out_cpt, out_tpr, top, ndx]
        }
        def next_xtcs = current_xtcs.mix(step_out.xtc)
        def next_sys = (c_idx == 1) ? step_out.sys_files : current_sys
        
        [next_state, next_xtcs, next_sys]
    }

    ch_all_xtcs = acc[1]
    ch_sys_files = acc[2]

    // 4. Group all resulting chunked XTC files by replicate ID
    ch_grouped_xtcs = ch_all_xtcs.groupTuple(by: 0, size: total_c)

    // 5. Concatenate raw chunks and expose for post-analysis
    TRJCAT(ch_grouped_xtcs)

    emit:
    raw_trajectory = TRJCAT.out.concat_xtc
    sys_files = ch_sys_files
}
