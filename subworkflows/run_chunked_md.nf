/*
 * subworkflows/run_chunked_md.nf
 * Replicated production MD where each 10 ns chunk is an independent Slurm job,
 * allowing faster cluster nodes to dynamically pick up subsequent chunks.
 */

include { GMX_EQUILIBRATION } from '../modules/gromacs/equilibration'
include { RUN_MD_INIT              } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_2  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_3  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_4  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_5  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_6  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_7  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_8  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_9  } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_10 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_11 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_12 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_13 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_14 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_15 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_16 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_17 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_18 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_19 } from '../modules/gromacs/run_chunk'
include { RUN_MD_CHUNK as CHUNK_20 } from '../modules/gromacs/run_chunk'
include { TRJCAT                   } from '../modules/gromacs/trjcat'
include { TRJCONV                  } from '../modules/gromacs/trjconv'

workflow RUN_CHUNKED_MD {
    take:
    ch_system        // channel: [ val(meta), path(gro), path(topol), path(ndx), path(posre_files) ]
    is_membrane      // val: true/false
    temperature      // val: 310.15
    total_ns         // val: 60
    chunk_ns         // val: 10
    replicates       // val: 3

    main:
    def num_chunks = Math.max(1, (total_ns / chunk_ns) as int)

    // 1. Equilibration on Slurm GPU
    GMX_EQUILIBRATION(
        ch_system,
        is_membrane,
        temperature
    )

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

    // 3. Chunk 1 (Initializes velocities & generates reference md_1.tpr)
    RUN_MD_INIT(ch_chunk1_in, is_membrane, chunk_ns, temperature)
    ch_all_xtcs = RUN_MD_INIT.out.xtc
    ch_ref_tpr  = RUN_MD_INIT.out.next_chunk_input.map { meta, gro, cpt, tpr -> [ meta, tpr ] }

    // Dynamically chained atomic chunks up to 200 ns
    ch_c2_in = RUN_MD_INIT.out.next_chunk_input.filter { num_chunks >= 2 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 2 ] }
    CHUNK_2(ch_c2_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_2.out.xtc)

    ch_c3_in = CHUNK_2.out.next_chunk_input.filter { num_chunks >= 3 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 3 ] }
    CHUNK_3(ch_c3_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_3.out.xtc)

    ch_c4_in = CHUNK_3.out.next_chunk_input.filter { num_chunks >= 4 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 4 ] }
    CHUNK_4(ch_c4_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_4.out.xtc)

    ch_c5_in = CHUNK_4.out.next_chunk_input.filter { num_chunks >= 5 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 5 ] }
    CHUNK_5(ch_c5_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_5.out.xtc)

    ch_c6_in = CHUNK_5.out.next_chunk_input.filter { num_chunks >= 6 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 6 ] }
    CHUNK_6(ch_c6_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_6.out.xtc)

    ch_c7_in = CHUNK_6.out.next_chunk_input.filter { num_chunks >= 7 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 7 ] }
    CHUNK_7(ch_c7_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_7.out.xtc)

    ch_c8_in = CHUNK_7.out.next_chunk_input.filter { num_chunks >= 8 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 8 ] }
    CHUNK_8(ch_c8_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_8.out.xtc)

    ch_c9_in = CHUNK_8.out.next_chunk_input.filter { num_chunks >= 9 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 9 ] }
    CHUNK_9(ch_c9_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_9.out.xtc)

    ch_c10_in = CHUNK_9.out.next_chunk_input.filter { num_chunks >= 10 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 10 ] }
    CHUNK_10(ch_c10_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_10.out.xtc)

    ch_c11_in = CHUNK_10.out.next_chunk_input.filter { num_chunks >= 11 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 11 ] }
    CHUNK_11(ch_c11_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_11.out.xtc)

    ch_c12_in = CHUNK_11.out.next_chunk_input.filter { num_chunks >= 12 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 12 ] }
    CHUNK_12(ch_c12_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_12.out.xtc)

    ch_c13_in = CHUNK_12.out.next_chunk_input.filter { num_chunks >= 13 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 13 ] }
    CHUNK_13(ch_c13_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_13.out.xtc)

    ch_c14_in = CHUNK_13.out.next_chunk_input.filter { num_chunks >= 14 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 14 ] }
    CHUNK_14(ch_c14_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_14.out.xtc)

    ch_c15_in = CHUNK_14.out.next_chunk_input.filter { num_chunks >= 15 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 15 ] }
    CHUNK_15(ch_c15_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_15.out.xtc)

    ch_c16_in = CHUNK_15.out.next_chunk_input.filter { num_chunks >= 16 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 16 ] }
    CHUNK_16(ch_c16_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_16.out.xtc)

    ch_c17_in = CHUNK_16.out.next_chunk_input.filter { num_chunks >= 17 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 17 ] }
    CHUNK_17(ch_c17_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_17.out.xtc)

    ch_c18_in = CHUNK_17.out.next_chunk_input.filter { num_chunks >= 18 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 18 ] }
    CHUNK_18(ch_c18_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_18.out.xtc)

    ch_c19_in = CHUNK_18.out.next_chunk_input.filter { num_chunks >= 19 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 19 ] }
    CHUNK_19(ch_c19_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_19.out.xtc)

    ch_c20_in = CHUNK_19.out.next_chunk_input.filter { num_chunks >= 20 }.map { meta, gro, cpt, tpr -> [ meta, cpt, tpr, 20 ] }
    CHUNK_20(ch_c20_in, chunk_ns)
    ch_all_xtcs = ch_all_xtcs.mix(CHUNK_20.out.xtc)

    // 4. Concatenate all chunk trajectories
    ch_to_concat = ch_all_xtcs.groupTuple(by: 0)
    TRJCAT(ch_to_concat)

    // 5. Clean trajectory (PBC correction and protein-ligand centering)
    ch_trjconv_in = TRJCAT.out.concat_xtc
        .join(ch_ref_tpr, by: 0)
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
