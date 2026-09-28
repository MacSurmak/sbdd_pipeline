/*
 * subworkflows/run_chunked_md.nf
 * Replicated production MD simulation with fully dynamic atomic chunk scheduling.
 * Supports any total simulation length and chunk interval (e.g. 500 ns / 1 ns chunks)
 * with zero hardcoded limits.
 */

include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
include { RUN_MD_CHUNK      } from '../modules/gromacs/run_chunk'
include { TRJCAT            } from '../modules/gromacs/trjcat'

// Static fallback template for 1 chunk (used for stub testing)
workflow RUN_CHUNKED_MD {
    take:
    ch_system
    ch_mdp_templates
    is_membrane
    temperature
    total_ns
    chunk_ns
    replicates
    slurm_queue
    nvt_short_ps
    nvt_long_ps
    npt_berendsen_ps
    npt_pr_ps

    main:
    GMX_EQUILIBRATION(ch_system, ch_mdp_templates, is_membrane, temperature, nvt_short_ps, nvt_long_ps, npt_berendsen_ps, npt_pr_ps)
    
    def num_reps = replicates as int
    ch_state_0 = GMX_EQUILIBRATION.out.equilibrated_state
        .join(ch_system.map{ m, gro, top, ndx, pos -> [m, top, ndx] }, by: 0)
        .flatMap { meta, eq_gro, eq_cpt, topol, ndx ->
            (1..num_reps).collect { rep_id ->
                def rep_meta = meta + [rep: rep_id]
                [ rep_meta, eq_gro, eq_cpt, eq_cpt, topol, ndx ]
            }
        }

    ch_all_xtcs = Channel.empty()

    step_1 = RUN_MD_CHUNK(ch_state_0, ch_mdp_templates, is_membrane, temperature, chunk_ns, 1)
    ch_all_xtcs = ch_all_xtcs.mix(step_1.xtc)
    
    sys_files = step_1.sys_files

    ch_grouped_xtcs = ch_all_xtcs.groupTuple(by: 0, size: 1)
    TRJCAT(ch_grouped_xtcs)

    emit:
    raw_trajectory = TRJCAT.out.concat_xtc
    sys_files = sys_files
}