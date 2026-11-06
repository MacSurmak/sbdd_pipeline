/*
 * subworkflows/dock_library.nf
 * Multi-Engine Docking and Pose Rescoring without embedded Prime MM-GBSA
 */

include { UNIDOCK_DOCK }            from '../modules/unidock/dock'
include { GNINA_RESCORE as GNINA_RESCORE_GLIDE }   from '../modules/gnina/rescore'
include { GNINA_RESCORE as GNINA_RESCORE_UNIDOCK } from '../modules/gnina/rescore'
include { GLIDE_GRIDGEN }           from '../modules/schrodinger/gridgen'
include { GLIDE_DOCK }              from '../modules/schrodinger/glide_dock'
include { ASSEMBLE_DOCKING_RECEPTOR; GENERATE_DOCKING_BOX; CALC_RMSD as CALC_RMSD_GLIDE; CALC_RMSD as CALC_RMSD_UNIDOCK; MERGE_UNIDOCK_DOCK as MERGE_UNIDOCK_CHUNKS; MERGE_GLIDE as MERGE_GLIDE_CHUNKS } from '../modules/docking/helpers'

workflow DOCK_LIBRARY {
    take:
    ch_receptor_pdb     // channel: [ val(meta), path(receptor_pdb) ]
    ch_ref_lig_sdf      // channel: [ val(meta), path(ref_lig_sdf) ]
    ch_chunks_sdf       // channel: [ val(meta), path(chunk_sdf), val(chunk_idx) ]
    ch_receptor_mae     // channel: [ val(meta), path(receptor_mae) ]
    ch_ref_lig_mae      // channel: [ val(meta), path(ref_lig_mae) ]
    ch_chunks_mae       // channel: [ val(meta), path(chunk_mae) ]
    ch_ion_sdfs
    ch_ion_maes
    ch_cof_sdfs
    ch_cof_maes
    dock_engine         // val: "both", "unidock", "glide"
    unidock_search_mode // val: "balance"
    unidock_scoring     // val: "vina"
    unidock_refine_step // val: 3
    unidock_num_modes   // val: 9
    unidock_energy_range // val: 5.0
    unidock_seed        // val: 42
    autobox_add         // val: 4.0
    inner_box           // val: "10, 10, 10"
    outer_box           // val: "26, 26, 26"
    glide_forcefield    // val: "OPLS4"
    glide_precision     // val: "SP"
    poses_per_lig       // val: 3
    gnina_enabled       // val: true/false
    gnina_rescore_mode  // val: "score_only" or "minimize"
    cnn_model           // val: ""
    glide_postdock      // val: true/false
    glide_recep_vscale  // val: 1.0
    glide_lig_vscale    // val: 1.0
    glide_recep_ccut    // val: 0.25
    glide_lig_ccut      // val: 0.15
    is_self             // val: boolean
    outdir              // val: canonical project outdir string
    dock_run_name       // val: unique DockingRun identifier

    main:
    ch_summary_files = Channel.empty()
    ch_primary_poses = Channel.empty()
    ch_receipts      = Channel.empty()
    def pub_subpath  = dock_run_name ? "docking/${dock_run_name}" : "docking"
    def do_gnina_rescore = (gnina_enabled.toString() == 'true')

    // 0. Generate PyMOL docking box script for every receptor pocket
    GENERATE_DOCKING_BOX(
        ch_ref_lig_sdf,
        inner_box,
        outer_box,
        autobox_add
    )

    // 1. ASSEMBLE HOLO RECEPTOR ONCE PER FRAME (APO + IONS + COFACTORS -> PDB, MAE, PDBQT)
    ch_assemble_in = ch_receptor_pdb
        .join(ch_receptor_mae, by: 0)
        .join(ch_ion_sdfs, by: 0)
        .join(ch_ion_maes, by: 0)
        .join(ch_cof_sdfs, by: 0)
        .join(ch_cof_maes, by: 0)

    ASSEMBLE_DOCKING_RECEPTOR(ch_assemble_in)

    ch_assembled_pdb   = ASSEMBLE_DOCKING_RECEPTOR.out.assembled_receptors.map { meta, pdb, mae, pdbqt -> [meta, pdb] }
    ch_assembled_mae   = ASSEMBLE_DOCKING_RECEPTOR.out.assembled_receptors.map { meta, pdb, mae, pdbqt -> [meta, mae] }
    ch_assembled_pdbqt = ASSEMBLE_DOCKING_RECEPTOR.out.assembled_receptors.map { meta, pdb, mae, pdbqt -> [meta, pdbqt] }
    ch_rec_rescore     = ch_assembled_pdb.join(ch_ref_lig_sdf, by: 0)

    // 2. GLIDE DOCKING (CPU) & GNINA 3D-CNN RESCORE
    if (dock_engine in ["glide", "both"]) {
        ch_grid_in = ch_assembled_mae.join(ch_ref_lig_mae, by: 0)
        GLIDE_GRIDGEN(
            ch_grid_in,
            inner_box,
            outer_box,
            glide_forcefield,
            glide_recep_vscale,
            glide_recep_ccut
        )

        ch_glide_in = GLIDE_GRIDGEN.out.grid_zip.combine(ch_chunks_mae)
            .filter { meta_g, grid_zip, meta_l, lib_chunk_mae ->
                if (!is_self) return true
                return meta_l.receptor == meta_g.id
            }
            .map { meta_g, grid_zip, meta_l, lib_chunk_mae ->
                def comb_meta = [
                    id: "${meta_g.id}_c${meta_l.chunk_idx}",
                    frame_id: meta_g.id,
                    frame: meta_g.frame,
                    chunk_idx: meta_l.chunk_idx,
                    project: meta_g.project
                ]
                [ comb_meta, grid_zip, lib_chunk_mae ]
            }
        GLIDE_DOCK(
            ch_glide_in,
            glide_precision,
            poses_per_lig,
            glide_forcefield,
            glide_postdock,
            glide_lig_vscale,
            glide_lig_ccut
        )

        ch_glide_pvs = GLIDE_DOCK.out.pv_maegz
            .map { comb_meta, pv -> [ [id: comb_meta.frame_id, frame: comb_meta.frame, project: comb_meta.project, run_name: dock_run_name], pv ] }
            .groupTuple(by: 0, sort: true)
        ch_glide_csvs = GLIDE_DOCK.out.scores_csv
            .map { comb_meta, csv -> [ [id: comb_meta.frame_id, frame: comb_meta.frame, project: comb_meta.project, run_name: dock_run_name], csv ] }
            .groupTuple(by: 0, sort: true)
        ch_glide_sdfs = GLIDE_DOCK.out.poses_sdf
            .map { comb_meta, sdf -> [ [id: comb_meta.frame_id, frame: comb_meta.frame, project: comb_meta.project, run_name: dock_run_name], sdf ] }
            .groupTuple(by: 0, sort: true)

        ch_to_merge_glide = ch_glide_pvs
            .join(ch_glide_csvs, by: 0)
            .join(ch_glide_sdfs, by: 0)

        MERGE_GLIDE_CHUNKS(ch_to_merge_glide)

        ch_summary_files = ch_summary_files
            .mix(MERGE_GLIDE_CHUNKS.out.scores_csv.map { meta, csv -> csv })
            .mix(MERGE_GLIDE_CHUNKS.out.poses_sdf.map { meta, sdf -> sdf })

        if (do_gnina_rescore) {
            ch_rescore_glide_in = MERGE_GLIDE_CHUNKS.out.poses_sdf
                .map { frame_meta, poses_sdf -> [ frame_meta.id, frame_meta, poses_sdf ] }
                .join(ch_rec_rescore.map { meta_rec, pdb, ref_sdf -> [ meta_rec.id, pdb, ref_sdf ] }, by: 0)
                .map { frame_id, frame_meta, poses_sdf, rec_pdb, ref_sdf ->
                    [ frame_meta, rec_pdb, ref_sdf, poses_sdf ]
                }
            GNINA_RESCORE_GLIDE(
                ch_rescore_glide_in,
                "glide",
                gnina_rescore_mode,
                autobox_add,
                cnn_model
            )
            ch_summary_files = ch_summary_files.mix(GNINA_RESCORE_GLIDE.out.rescored_sdf.map { meta, sdf -> sdf })
            ch_glide_final_poses = GNINA_RESCORE_GLIDE.out.rescored_sdf
        } else {
            ch_glide_final_poses = MERGE_GLIDE_CHUNKS.out.poses_sdf
        }

        if (is_self) {
            ch_glide_rmsd_in = ch_glide_final_poses
                .map { meta, docked_sdf -> [ meta.id, meta, docked_sdf ] }
                .join(ch_ref_lig_sdf.map { meta, ref_sdf -> [ meta.id, ref_sdf ] }, by: 0)
                .map { frame_id, meta, docked_sdf, ref_sdf -> [ meta, docked_sdf, ref_sdf ] }

            CALC_RMSD_GLIDE(ch_glide_rmsd_in, 'glide')
            ch_summary_files = ch_summary_files.mix(CALC_RMSD_GLIDE.out.rmsd_csv.map { meta, csv -> csv })
        }

        ch_primary_poses = ch_glide_final_poses
        ch_receipts = ch_receipts.mix(
            ch_glide_final_poses.map { meta, final_sdf ->
                "${meta.id}\tglide\t${outdir}/${pub_subpath}/${meta.id}/${final_sdf.name}"
            }
        )
    }

    // 3. UNI-DOCK GPU DOCKING & GNINA 3D-CNN RESCORE
    if (dock_engine in ["unidock", "both"]) {
        ch_rec_unidock = ch_assembled_pdbqt.join(ch_ref_lig_sdf, by: 0)
        ch_unidock_in = ch_chunks_sdf.combine(ch_rec_unidock)
            .filter { meta_chunk, chunk_sdf, chunk_idx, meta_rec, rec_pdbqt, ref_sdf ->
                if (!is_self) return true
                return meta_chunk.receptor == meta_rec.id
            }
            .map { meta_chunk, chunk_sdf, chunk_idx, meta_rec, rec_pdbqt, ref_sdf ->
                def comb_meta = [
                    id: "${meta_rec.id}_c${chunk_idx}",
                    frame_id: meta_rec.id,
                    frame: meta_rec.frame,
                    project: meta_rec.project
                ]
                [ comb_meta, rec_pdbqt, ref_sdf, chunk_sdf, chunk_idx ]
            }

        UNIDOCK_DOCK(
            ch_unidock_in,
            unidock_search_mode,
            unidock_scoring,
            unidock_refine_step,
            unidock_num_modes,
            unidock_energy_range,
            autobox_add,
            unidock_seed
        )

        ch_to_merge_unidock = UNIDOCK_DOCK.out.docked_sdf
            .map { meta, sdf ->
                def frame_meta = [id: meta.frame_id, frame: meta.frame, project: meta.project, run_name: dock_run_name]
                [ frame_meta, sdf ]
            }
            .groupTuple(by: 0, sort: true)

        MERGE_UNIDOCK_CHUNKS(ch_to_merge_unidock)

        ch_summary_files = ch_summary_files
            .mix(MERGE_UNIDOCK_CHUNKS.out.merged_sdf.map { meta, sdf -> sdf })

        if (do_gnina_rescore) {
            ch_rescore_unidock_in = MERGE_UNIDOCK_CHUNKS.out.merged_sdf
                .map { frame_meta, poses_sdf -> [ frame_meta.id, frame_meta, poses_sdf ] }
                .join(ch_rec_rescore.map { meta_rec, pdb, ref_sdf -> [ meta_rec.id, pdb, ref_sdf ] }, by: 0)
                .map { frame_id, frame_meta, poses_sdf, rec_pdb, ref_sdf ->
                    [ frame_meta, rec_pdb, ref_sdf, poses_sdf ]
                }
            GNINA_RESCORE_UNIDOCK(
                ch_rescore_unidock_in,
                "unidock",
                gnina_rescore_mode,
                autobox_add,
                cnn_model
            )
            ch_summary_files = ch_summary_files.mix(GNINA_RESCORE_UNIDOCK.out.rescored_sdf.map { meta, sdf -> sdf })
            ch_unidock_final_poses = GNINA_RESCORE_UNIDOCK.out.rescored_sdf
        } else {
            ch_unidock_final_poses = MERGE_UNIDOCK_CHUNKS.out.merged_sdf
        }

        if (is_self) {
            ch_unidock_rmsd_in = ch_unidock_final_poses
                .map { meta, docked_sdf -> [ meta.id, meta, docked_sdf ] }
                .join(ch_ref_lig_sdf.map { meta, ref_sdf -> [ meta.id, ref_sdf ] }, by: 0)
                .map { frame_id, meta, docked_sdf, ref_sdf -> [ meta, docked_sdf, ref_sdf ] }

            CALC_RMSD_UNIDOCK(ch_unidock_rmsd_in, 'unidock')
            ch_summary_files = ch_summary_files.mix(CALC_RMSD_UNIDOCK.out.rmsd_csv.map { meta, csv -> csv })
        }

        if (dock_engine == "unidock") {
            ch_primary_poses = ch_unidock_final_poses
        }

        ch_receipts = ch_receipts.mix(
            ch_unidock_final_poses.map { meta, final_sdf ->
                "${meta.id}\tunidock\t${outdir}/${pub_subpath}/${meta.id}/${final_sdf.name}"
            }
        )
    }

    ch_all_summary_files = ch_summary_files.collect()

    emit:
    docked_poses    = ch_primary_poses
    summary_files   = ch_all_summary_files
    receipt_entries = ch_receipts
}
