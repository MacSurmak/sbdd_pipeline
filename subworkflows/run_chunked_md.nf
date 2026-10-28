/*
 * =========================================================================================
 * ARCHITECTURAL DEFECT & TECHNICAL DEBT NOTICE: STATIC CHUNK UNROLLING LIMIT
 * =========================================================================================
 * ISSUE:
 *   This subworkflow statically instantiates exactly 10 sequential MD chunk execution steps
 *   (RUN_MD_CHUNK_1 through RUN_MD_CHUNK_10). Nextflow DSL2 lacks native dynamic cyclic
 *   recursion without pre-compiling explicit process nodes into the DAG.
 *
 * CONSEQUENCE:
 *   If 'assets/configs/md.yaml' specifies total_time_ns > 10 * chunk_size_ns (e.g. 200 ns total
 *   with 10 ns chunks = 20 chunks), any chunks beyond #10 are silently omitted, or the simulation
 *   stops at 100 ns.
 *
 * REQUIRED FUNDAMENTAL RESOLUTION:
 *   Replace this statically unrolled subworkflow with an iterative Nextflow orchestrator pattern
 *   (e.g., dynamic external submission via Nextflow feedback channel, an orchestrating CLI loop
 *   in 'sbdd md' that submits batches of chunks, or a recursive DAG generator that dynamically
 *   templates the exact number of required chunk processes based on md.yaml before invocation).
 * =========================================================================================
 */
/* Dynamically generated unrolled DAG for chunked MD */
include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
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
include { RUN_MD_CHUNK as RUN_MD_CHUNK_11 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_12 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_13 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_14 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_15 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_16 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_17 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_18 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_19 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as RUN_MD_CHUNK_20 } from '../modules/gromacs/run_chunk'

workflow RUN_CHUNKED_MD {
    take:
    ch_system
    ch_md_protocol
    total_ns
    chunk_ns
    replicates
    minim_flags
    eq_flags
    prod_flags

    main:
    GMX_EQUILIBRATION(ch_system, ch_md_protocol, minim_flags, eq_flags)

    def num_reps = replicates as int
    def num_chunks = Math.round((total_ns as double) / (chunk_ns as double)) as int
    if (num_chunks < 1 || num_chunks > 20) {
        error "FATAL: Calculated MD chunk count (${num_chunks}) is outside supported range (1..20)!"
    }

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

    step_1 = RUN_MD_CHUNK_1(ch_state_0, ch_md_protocol, chunk_ns, 1, prod_flags)
    ch_curr = step_1.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
    ch_all_xtcs = ch_all_xtcs.mix(step_1.xtc)

    if (num_chunks >= 2) {
        step_2 = RUN_MD_CHUNK_2(ch_curr, ch_md_protocol, chunk_ns, 2, prod_flags)
        ch_curr = step_2.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_2.xtc)
    }
    if (num_chunks >= 3) {
        step_3 = RUN_MD_CHUNK_3(ch_curr, ch_md_protocol, chunk_ns, 3, prod_flags)
        ch_curr = step_3.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_3.xtc)
    }
    if (num_chunks >= 4) {
        step_4 = RUN_MD_CHUNK_4(ch_curr, ch_md_protocol, chunk_ns, 4, prod_flags)
        ch_curr = step_4.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_4.xtc)
    }
    if (num_chunks >= 5) {
        step_5 = RUN_MD_CHUNK_5(ch_curr, ch_md_protocol, chunk_ns, 5, prod_flags)
        ch_curr = step_5.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_5.xtc)
    }
    if (num_chunks >= 6) {
        step_6 = RUN_MD_CHUNK_6(ch_curr, ch_md_protocol, chunk_ns, 6, prod_flags)
        ch_curr = step_6.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_6.xtc)
    }
    if (num_chunks >= 7) {
        step_7 = RUN_MD_CHUNK_7(ch_curr, ch_md_protocol, chunk_ns, 7, prod_flags)
        ch_curr = step_7.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_7.xtc)
    }
    if (num_chunks >= 8) {
        step_8 = RUN_MD_CHUNK_8(ch_curr, ch_md_protocol, chunk_ns, 8, prod_flags)
        ch_curr = step_8.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_8.xtc)
    }
    if (num_chunks >= 9) {
        step_9 = RUN_MD_CHUNK_9(ch_curr, ch_md_protocol, chunk_ns, 9, prod_flags)
        ch_curr = step_9.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_9.xtc)
    }
    if (num_chunks >= 10) {
        step_10 = RUN_MD_CHUNK_10(ch_curr, ch_md_protocol, chunk_ns, 10, prod_flags)
        ch_curr = step_10.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_10.xtc)
    }
    if (num_chunks >= 11) {
        step_11 = RUN_MD_CHUNK_11(ch_curr, ch_md_protocol, chunk_ns, 11, prod_flags)
        ch_curr = step_11.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_11.xtc)
    }
    if (num_chunks >= 12) {
        step_12 = RUN_MD_CHUNK_12(ch_curr, ch_md_protocol, chunk_ns, 12, prod_flags)
        ch_curr = step_12.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_12.xtc)
    }
    if (num_chunks >= 13) {
        step_13 = RUN_MD_CHUNK_13(ch_curr, ch_md_protocol, chunk_ns, 13, prod_flags)
        ch_curr = step_13.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_13.xtc)
    }
    if (num_chunks >= 14) {
        step_14 = RUN_MD_CHUNK_14(ch_curr, ch_md_protocol, chunk_ns, 14, prod_flags)
        ch_curr = step_14.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_14.xtc)
    }
    if (num_chunks >= 15) {
        step_15 = RUN_MD_CHUNK_15(ch_curr, ch_md_protocol, chunk_ns, 15, prod_flags)
        ch_curr = step_15.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_15.xtc)
    }
    if (num_chunks >= 16) {
        step_16 = RUN_MD_CHUNK_16(ch_curr, ch_md_protocol, chunk_ns, 16, prod_flags)
        ch_curr = step_16.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_16.xtc)
    }
    if (num_chunks >= 17) {
        step_17 = RUN_MD_CHUNK_17(ch_curr, ch_md_protocol, chunk_ns, 17, prod_flags)
        ch_curr = step_17.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_17.xtc)
    }
    if (num_chunks >= 18) {
        step_18 = RUN_MD_CHUNK_18(ch_curr, ch_md_protocol, chunk_ns, 18, prod_flags)
        ch_curr = step_18.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_18.xtc)
    }
    if (num_chunks >= 19) {
        step_19 = RUN_MD_CHUNK_19(ch_curr, ch_md_protocol, chunk_ns, 19, prod_flags)
        ch_curr = step_19.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_19.xtc)
    }
    if (num_chunks >= 20) {
        step_20 = RUN_MD_CHUNK_20(ch_curr, ch_md_protocol, chunk_ns, 20, prod_flags)
        ch_curr = step_20.state.map { meta, out_gro, out_cpt, out_tpr, top, ndx, pos -> [meta, out_gro, out_cpt, out_tpr, top, ndx, pos] }
        ch_all_xtcs = ch_all_xtcs.mix(step_20.xtc)
    }

    sys_files = step_1.sys_files

    emit:
    sys_files = sys_files
}
