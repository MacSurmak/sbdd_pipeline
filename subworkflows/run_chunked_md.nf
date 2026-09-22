/*
 * subworkflows/run_chunked_md.nf
 * Replicated production MD simulation with fully dynamic atomic chunk scheduling.
 * Supports any total simulation length and chunk interval (e.g. 500 ns / 1 ns chunks)
 * with zero hardcoded limits.
 */

include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
include { RUN_PRODUCTION_MD  } from '../modules/gromacs/run_chunk'

workflow RUN_CHUNKED_MD {
    take:
    ch_system        // channel: [ val(meta), path(gro), path(topol), path(ndx), path(posre_files) ]
    is_membrane      // val: true/false
    temperature      // val: 310.15
    total_ns         // val: 60
    chunk_ns         // val: 10
    replicates       // val: 3
    slurm_queue      // val: "main"

    main:
    // 1. Equilibration on Slurm GPU with dynamic step counts
    GMX_EQUILIBRATION(
        ch_system,
        is_membrane,
        temperature
    )

    ch_sys_info = ch_system.map { meta, gro, topol, ndx, posre -> [ meta.id, topol, ndx ] }

    // 2. Expand into independent replicates for production MD
    def num_reps = (replicates as int)
    ch_prod_in = GMX_EQUILIBRATION.out.equilibrated_state
        .map { meta, eq_gro, eq_cpt -> [ meta.id, meta, eq_gro, eq_cpt ] }
        .combine(ch_sys_info, by: 0)
        .flatMap { target_id, meta, eq_gro, eq_cpt, topol, ndx ->
            (1..num_reps).collect { rep_id ->
                def rep_meta = meta + [rep: rep_id]
                [ rep_meta, eq_gro, eq_cpt, topol, ndx ]
            }
        }

    // 3. Dynamic Production MD
    RUN_PRODUCTION_MD(
        ch_prod_in,
        is_membrane,
        temperature,
        total_ns,
        chunk_ns,
        slurm_queue
    )

    emit:
    clean_trajectory = RUN_PRODUCTION_MD.out.clean_traj
}
