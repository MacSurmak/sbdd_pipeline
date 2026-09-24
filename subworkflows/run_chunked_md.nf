/*
 * subworkflows/run_chunked_md.nf
 * Replicated production MD simulation with fully dynamic atomic chunk scheduling.
 * Supports any total simulation length and chunk interval (e.g. 500 ns / 1 ns chunks)
 * with zero hardcoded limits.
 */

include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
include { RUN_MD_CHUNK      } from '../modules/gromacs/run_chunk'
include { TRJCAT            } from '../modules/gromacs/trjcat'

// Pure Nextflow Recursive Loop for Slurm Chunk Offloading
workflow RECURSIVE_MD {
    take:
    state_in // [ meta, chunk_idx, gro, cpt, tpr, top, ndx ]
    mdp_templates
    is_membrane
    temperature
    chunk_ns

    main:
    RUN_MD_CHUNK(state_in, mdp_templates, is_membrane, temperature, chunk_ns)

    RUN_MD_CHUNK.out.state.branch { meta, c_idx, gro, cpt, tpr, top, ndx ->
        done: c_idx >= meta.total_chunks
        loop: c_idx < meta.total_chunks
    }.set { chunk_router }

    next_state = chunk_router.loop.map { meta, c_idx, gro, cpt, tpr, top, ndx ->
        [meta, c_idx + 1, gro, cpt, tpr, top, ndx]
    }

    RECURSIVE_MD(next_state, mdp_templates, is_membrane, temperature, chunk_ns)

    emit:
    xtcs = RUN_MD_CHUNK.out.xtc
    sys_files = RUN_MD_CHUNK.out.sys_files
}

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

    // 2. Expand into independent replicates for dynamic recursive production MD
    ch_start_state = GMX_EQUILIBRATION.out.equilibrated_state
        .join(ch_system.map{ m, gro, top, ndx, pos -> [m.id, top, ndx] }, by: 0)
        .flatMap { target_id, meta, eq_gro, eq_cpt, topol, ndx ->
            (1..num_reps).collect { rep_id ->
                def rep_meta = meta + [rep: rep_id, total_chunks: total_c]
                // Send eq_cpt in place of prev_tpr as a dummy for the very first chunk
                [ rep_meta, 1, eq_gro, eq_cpt, eq_cpt, topol, ndx ]
            }
        }

    // 3. Initiate self-balancing recursive loop
    RECURSIVE_MD(ch_start_state, ch_mdp_templates, is_membrane, temperature, chunk_ns)

    // 4. Group all resulting chunked XTC files by replicate ID
    ch_xtcs_grouped = RECURSIVE_MD.out.xtcs.groupTuple(by: 0)

    // 5. Concatenate raw chunks and expose for post-analysis
    TRJCAT(ch_xtcs_grouped)

    emit:
    raw_trajectory = TRJCAT.out.concat_xtc
    sys_files = RECURSIVE_MD.out.sys_files
}
