/* Dynamically generated unrolled DAG for chunked MD */
include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
include { TRJCAT            } from '../modules/gromacs/trjcat'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_1 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_2 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_3 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_4 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_5 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_6 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_7 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_8 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_9 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_10 } from '../modules/gromacs/run_chunk'

workflow RUN_CHUNKED_MD {
    take:
    ch_system
    ch_md_protocol
    total_ns
    chunk_ns
    replicates

    main:
    GMX_EQUILIBRATION(ch_system, ch_md_protocol)
    
    def num_reps = replicates as int
    ch_state_0 = GMX_EQUILIBRATION.out.equilibrated_state
        .join(ch_system.map{ m, gro, top, ndx, pos -> [m, top, ndx, pos] }, by: 0)
        .flatMap { meta, eq_gro, eq_cpt, eq_tpr, topol, ndx, posre_files ->
            def rep_ids = meta.replicates ?: (1..num_reps).toList()
            rep_ids.collect { rep_id ->
                def rep_meta = meta + [rep: rep_id]
                [ rep_meta, eq_gro, eq_cpt, eq_tpr, topol, ndx, posre_files ]
            }
        }

    ch_all_xtcs = Channel.empty()

    step_1 = RUN_MD_CHUNK_1(ch_state_0, ch_md_protocol, chunk_ns, 1)
    ch_state_1 = step_1.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_1.xtc)

    step_2 = RUN_MD_CHUNK_2(ch_state_1, ch_md_protocol, chunk_ns, 2)
    ch_state_2 = step_2.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_2.xtc)

    step_3 = RUN_MD_CHUNK_3(ch_state_2, ch_md_protocol, chunk_ns, 3)
    ch_state_3 = step_3.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_3.xtc)

    step_4 = RUN_MD_CHUNK_4(ch_state_3, ch_md_protocol, chunk_ns, 4)
    ch_state_4 = step_4.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_4.xtc)

    step_5 = RUN_MD_CHUNK_5(ch_state_4, ch_md_protocol, chunk_ns, 5)
    ch_state_5 = step_5.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_5.xtc)

    step_6 = RUN_MD_CHUNK_6(ch_state_5, ch_md_protocol, chunk_ns, 6)
    ch_state_6 = step_6.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_6.xtc)

    step_7 = RUN_MD_CHUNK_7(ch_state_6, ch_md_protocol, chunk_ns, 7)
    ch_state_7 = step_7.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_7.xtc)

    step_8 = RUN_MD_CHUNK_8(ch_state_7, ch_md_protocol, chunk_ns, 8)
    ch_state_8 = step_8.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_8.xtc)

    step_9 = RUN_MD_CHUNK_9(ch_state_8, ch_md_protocol, chunk_ns, 9)
    ch_state_9 = step_9.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_9.xtc)

    step_10 = RUN_MD_CHUNK_10(ch_state_9, ch_md_protocol, chunk_ns, 10)
    ch_state_10 = step_10.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_10.xtc)

    sys_files = step_1.sys_files

    ch_grouped_xtcs = ch_all_xtcs.groupTuple(by: 0, size: 10, sort: true)
    TRJCAT(ch_grouped_xtcs)

    emit:
    raw_trajectory = TRJCAT.out.concat_xtc
    sys_files = sys_files
}
