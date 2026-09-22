/*
 * subworkflows/run_chunked_md.nf
 * Replicated production MD where each 10 ns chunk is an independent Slurm job,
 * allowing faster cluster nodes to dynamically pick up subsequent chunks.
 */

include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
include { RUN_MD_INIT       } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_2 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_3 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_4 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_5 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_6 } from '../modules/gromacs/run_chunk'
include { TRJCAT            } from '../modules/gromacs/trjcat'
include { TRJCONV           } from '../modules/gromacs/trjconv'

workflow RUN_CHUNKED_MD {
    take:
    ch_system        // channel: [ val(meta), path(gro), path(topol), path(ndx), path(posre_files) ]
    is_membrane      // val: true/false
    temperature      // val: 310.15
    total_ns         // val: 60
    chunk_ns         // val: 10
    replicates       // val: 3

    main:
    // 1. Equilibration on Slurm GPU
    GMX_EQUILIBRATION(
        ch_system,
        is_membrane,
        temperature
    )

    // System topology/index lookup keyed by target_id
    ch_sys_info = ch_system.map { meta, gro, topol, ndx, posre -> [ meta.id, topol, ndx ] }

    // 2. Expand into independent replicates for production MD
    def num_reps = (replicates as int)
    ch_chunk1_in = GMX_EQUILIBRATION.out.equilibrated_state
        .map { meta, eq_gro, eq_cpt -> [ meta.id, meta, eq_gro, eq_cpt ] }
        .combine(ch_sys_info, by: 0)
        .flatMap { target_id, meta, eq_gro, eq_cpt, topol, ndx ->
            (1..num_reps).collect { rep_id ->
                def rep_meta = meta + [rep: rep_id]
                [ rep_meta, eq_gro, eq_cpt, topol, ndx ]
            }
        }

    // 3. Chunk 1 (0 -> 10 ns) - Initializes velocities & generates md_1.tpr
    RUN_MD_INIT(ch_chunk1_in, is_membrane, chunk_ns, temperature)
    ch_all_xtcs = RUN_MD_INIT.out.xtc
    ch_last_tpr = RUN_MD_INIT.out.next_chunk_input.map { meta, gro, cpt, tpr -> [ meta, tpr ] }

    // Chunk 2 (10 -> 20 ns)
    ch_chunk2_in = RUN_MD_INIT.out.next_chunk_input
        .map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 2 ] }
    CHUNK_2(ch_chunk2_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_2.out.xtc)
    ch_last_tpr = CHUNK_2.out.next_chunk_input.map { meta, gro, cpt, tpr -> [ meta, tpr ] }

    // Chunk 3 (20 -> 30 ns)
    ch_chunk3_in = CHUNK_2.out.next_chunk_input
        .map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 3 ] }
    CHUNK_3(ch_chunk3_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_3.out.xtc)
    ch_last_tpr = CHUNK_3.out.next_chunk_input.map { meta, gro, cpt, tpr -> [ meta, tpr ] }

    // Chunk 4 (30 -> 40 ns)
    ch_chunk4_in = CHUNK_3.out.next_chunk_input
        .map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 4 ] }
    CHUNK_4(ch_chunk4_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_4.out.xtc)
    ch_last_tpr = CHUNK_4.out.next_chunk_input.map { meta, gro, cpt, tpr -> [ meta, tpr ] }

    // Chunk 5 (40 -> 50 ns)
    ch_chunk5_in = CHUNK_4.out.next_chunk_input
        .map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 5 ] }
    CHUNK_5(ch_chunk5_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_5.out.xtc)
    ch_last_tpr = CHUNK_5.out.next_chunk_input.map { meta, gro, cpt, tpr -> [ meta, tpr ] }

    // Chunk 6 (50 -> 60 ns)
    ch_chunk6_in = CHUNK_5.out.next_chunk_input
        .map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 6 ] }
    CHUNK_6(ch_chunk6_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_6.out.xtc)
    ch_last_tpr = CHUNK_6.out.next_chunk_input.map { meta, gro, cpt, tpr -> [ meta, tpr ] }

    // 4. Concatenate all chunk xtcs for each replicate
    ch_to_concat = ch_all_xtcs.groupTuple(by: 0)
    TRJCAT(ch_to_concat)

    // 5. Clean trajectory (PBC fix + centering)
    ch_trjconv_in = TRJCAT.out.concat_xtc
        .join(ch_last_tpr, by: 0)
        .map { meta, xtc, tpr -> [ meta.id, meta, xtc, tpr ] }
        .combine(ch_system.map { meta, gro, topol, ndx, posre -> [ meta.id, ndx ] }, by: 0)
        .map { target_id, meta, xtc, tpr, ndx -> [ meta, xtc, tpr, ndx ] }

    TRJCONV(ch_trjconv_in)

    // Join with topol and ndx for downstream analysis
    ch_ready_analysis = TRJCONV.out.clean_traj
        .map { meta, clean_xtc, ref_pdb -> [ meta.id, meta, clean_xtc, ref_pdb ] }
        .combine(ch_system.map { meta, gro, topol, ndx, posre -> [ meta.id, topol, ndx ] }, by: 0)
        .map { target_id, meta, clean_xtc, ref_pdb, topol, ndx ->
            [ meta, clean_xtc, ref_pdb, topol, ndx ]
        }

    emit:
    clean_trajectory = ch_ready_analysis
}
