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

    ch_all_xtcs = Channel.empty()
    ch_sys_files = Channel.empty()

    // 3. Initiate self-balancing unrolled DAG loop
    for (int i = 1; i <= total_c; i++) {
        def step_out = RUN_MD_CHUNK(ch_state, ch_mdp_templates, is_membrane, temperature, chunk_ns, i)
        ch_state = step_out.state
        ch_all_xtcs = ch_all_xtcs.mix(step_out.xtc)
        if (i == 1) {
            ch_sys_files = step_out.sys_files
        }
    }

    // 4. Group all resulting chunked XTC files by replicate ID (strictly waiting for all chunks)
    ch_grouped_xtcs = ch_all_xtcs.groupTuple(by: 0, size: total_c)

    // 5. Concatenate raw chunks and expose for post-analysis
    TRJCAT(ch_grouped_xtcs)

    emit:
    raw_trajectory = TRJCAT.out.concat_xtc
    sys_files = ch_sys_files
}
