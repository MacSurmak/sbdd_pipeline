/*
 * subworkflows/dock_library.nf
 * Multi-Engine Docking and Pose Rescoring without embedded Prime MM-GBSA
 */

include { GNINA_DOCK }              from '../modules/gnina/dock'
include { GNINA_RESCORE }           from '../modules/gnina/rescore'
include { GLIDE_GRIDGEN }           from '../modules/schrodinger/gridgen'
include { GLIDE_DOCK }              from '../modules/schrodinger/glide_dock'
include { GENERATE_DOCKING_BOX; CALC_RMSD as CALC_RMSD_GLIDE; CALC_RMSD as CALC_RMSD_GNINA; MERGE_SDF as MERGE_DOCK_CHUNKS; MERGE_SDF as MERGE_RESC_CHUNKS; MERGE_GLIDE as MERGE_GLIDE_CHUNKS } from '../modules/docking/helpers'

workflow DOCK_LIBRARY {
    take:
    ch_receptor_pdb     // channel: [ val(meta), path(receptor_pdb) ]
    ch_ref_lig_sdf      // channel: [ val(meta), path(ref_lig_sdf) ]
    ch_chunks_sdf       // channel: [ val(meta), path(chunk_sdf), val(chunk_idx) ]
    ch_receptor_mae     // channel: [ val(meta), path(receptor_mae) ]
    ch_ref_lig_mae      // channel: [ val(meta), path(ref_lig_mae) ]
    ch_chunks_mae       // channel: [ val(meta), path(chunk_mae) ]
    dock_engine         // val: "both", "gnina", "glide"
    exhaustiveness      // val: 16
    cnn_scoring         // val: "rescore"
    num_modes           // val: 3
    autobox_add         // val: 4.0
    is_self             // val: boolean

    main:
    ch_summary_files = Channel.empty()
    ch_primary_poses = Channel.empty()

    // 0. Generate PyMOL docking box script for every receptor pocket
    GENERATE_DOCKING_BOX(
        ch_ref_lig_sdf,
        params.inner_box,
        params.outer_box,
        autobox_add
    )
    ch_summary_files = ch_summary_files.mix(GENERATE_DOCKING_BOX.out.box_pml.map { meta, pml -> pml })

    // 1. GLIDE DOCKING & GNINA RESCORE
    if (dock_engine in ["glide", "both"]) {
        def inner_box = params.inner_box
        def outer_box = params.outer_box
        def ff_dock   = params.glide_forcefield
        def prec      = params.glide_precision
        def poses_n   = params.poses_per_lig
        def cnn_resc  = params.cnn_model
        def r_vscale  = params.glide_recep_vscale ?: 0.85
        def l_vscale  = params.glide_lig_vscale ?: 0.85

        ch_grid_in = ch_receptor_mae.join(ch_ref_lig_mae, by: 0)
        GLIDE_GRIDGEN(ch_grid_in, inner_box, outer_box, ff_dock, r_vscale)

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
        GLIDE_DOCK(ch_glide_in, prec, poses_n, ff_dock, l_vscale)

        ch_rec_rescore = ch_receptor_pdb.join(ch_ref_lig_sdf, by: 0)
        ch_rescore_in = GLIDE_DOCK.out.poses_sdf
            .map { comb_meta, sdf -> [ comb_meta.frame_id, comb_meta, sdf ] }
            .combine(ch_rec_rescore.map { meta_rec, pdb, ref_sdf -> [ meta_rec.id, pdb, ref_sdf ] }, by: 0)
            .map { frame_id, comb_meta, poses_sdf, rec_pdb, ref_sdf ->
                [ comb_meta, rec_pdb, ref_sdf, poses_sdf ]
            }
        GNINA_RESCORE(ch_rescore_in, cnn_resc)

        // Consolidate Glide chunks per frame
        ch_to_merge_glide = GLIDE_DOCK.out.pv_maegz
            .join(GLIDE_DOCK.out.scores_csv, by: 0)
            .join(GLIDE_DOCK.out.poses_sdf, by: 0)
            .map { comb_meta, pv, csv, sdf ->
                def frame_meta = [
                    id: comb_meta.frame_id,
                    frame: comb_meta.frame,
                    project: comb_meta.project
                ]
                [ frame_meta, pv, csv, sdf ]
            }
            .groupTuple(by: 0, sort: true)

        MERGE_GLIDE_CHUNKS(ch_to_merge_glide)

        // Consolidate Rescore chunks per frame
        ch_to_merge_resc = GNINA_RESCORE.out.rescored_sdf
            .map { comb_meta, sdf ->
                def frame_meta = [
                    id: comb_meta.frame_id,
                    frame: comb_meta.frame,
                    project: comb_meta.project
                ]
                [ frame_meta, sdf ]
            }
            .groupTuple(by: 0, sort: true)

        MERGE_RESC_CHUNKS(ch_to_merge_resc)

        // In-situ RMSD calculation for Glide poses
        ch_glide_rmsd_in = MERGE_GLIDE_CHUNKS.out.poses_sdf
            .map { meta, docked_sdf -> [ meta.id, meta, docked_sdf ] }
            .join(ch_ref_lig_sdf.map { meta, ref_sdf -> [ meta.id, ref_sdf ] }, by: 0)
            .map { frame_id, meta, docked_sdf, ref_sdf -> [ meta, docked_sdf, ref_sdf ] }

        CALC_RMSD_GLIDE(ch_glide_rmsd_in, 'glide')

        ch_summary_files = ch_summary_files
            .mix(MERGE_GLIDE_CHUNKS.out.scores_csv.map { meta, csv -> csv })
            .mix(MERGE_GLIDE_CHUNKS.out.poses_sdf.map { meta, sdf -> sdf })
            .mix(MERGE_RESC_CHUNKS.out.merged_sdf.map { meta, sdf -> sdf })
            .mix(CALC_RMSD_GLIDE.out.rmsd_csv.map { meta, csv -> csv })
            .mix(ch_ref_lig_sdf.map { meta, sdf -> sdf })

        ch_primary_poses = MERGE_GLIDE_CHUNKS.out.poses_sdf
    }

    // 2. GNINA DOCKING
    if (dock_engine in ["gnina", "both"]) {
        ch_rec_comb = ch_receptor_pdb.join(ch_ref_lig_sdf, by: 0)
        ch_dock_in = ch_chunks_sdf.combine(ch_rec_comb)
            .filter { meta_chunk, chunk_sdf, chunk_idx, meta_rec, rec_pdb, ref_sdf ->
                if (!is_self) return true
                return meta_chunk.receptor == meta_rec.id
            }
            .map { meta_chunk, chunk_sdf, chunk_idx, meta_rec, rec_pdb, ref_sdf ->
                def comb_meta = [
                    id: "${meta_rec.id}_c${chunk_idx}",
                    frame_id: meta_rec.id,
                    frame: meta_rec.frame,
                    project: meta_rec.project
                ]
                [ comb_meta, rec_pdb, ref_sdf, chunk_sdf, chunk_idx ]
            }

        GNINA_DOCK(
            ch_dock_in,
            exhaustiveness,
            cnn_scoring,
            num_modes,
            autobox_add
        )

        ch_to_merge = GNINA_DOCK.out.docked_sdf
            .map { meta, sdf ->
                def frame_meta = [id: meta.frame_id, frame: meta.frame, project: meta.project]
                [ frame_meta, sdf ]
            }
            .groupTuple(by: 0, sort: true)

        MERGE_DOCK_CHUNKS(ch_to_merge)

        // In-situ RMSD calculation for GNINA poses
        ch_gnina_rmsd_in = MERGE_DOCK_CHUNKS.out.merged_sdf
            .map { meta, docked_sdf -> [ meta.id, meta, docked_sdf ] }
            .join(ch_ref_lig_sdf.map { meta, ref_sdf -> [ meta.id, ref_sdf ] }, by: 0)
            .map { frame_id, meta, docked_sdf, ref_sdf -> [ meta, docked_sdf, ref_sdf ] }

        CALC_RMSD_GNINA(ch_gnina_rmsd_in, 'gnina')

        ch_summary_files = ch_summary_files
            .mix(MERGE_DOCK_CHUNKS.out.merged_sdf.map { meta, sdf -> sdf })
            .mix(CALC_RMSD_GNINA.out.rmsd_csv.map { meta, csv -> csv })

        if (dock_engine == "gnina") {
            ch_primary_poses = MERGE_DOCK_CHUNKS.out.merged_sdf
        }
    }

    ch_all_summary_files = ch_summary_files.collect()

    emit:
    docked_poses   = ch_primary_poses
    summary_files  = ch_all_summary_files
}
