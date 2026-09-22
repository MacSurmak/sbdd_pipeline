/*
 * subworkflows/dock_library.nf
 * Multi-Engine Docking, Rescoring, Funnel Triage, and Prime MM-GBSA
 */

include { GNINA_DOCK }              from '../modules/gnina/dock'
include { GNINA_RESCORE }           from '../modules/gnina/rescore'
include { GLIDE_GRIDGEN }           from '../modules/schrodinger/gridgen'
include { GLIDE_DOCK }              from '../modules/schrodinger/glide_dock'
include { BUILD_FUNNEL_WHITELIST }  from '../modules/reporting/build_funnel'
include { FILTER_PV }               from '../modules/schrodinger/filter_pv'
include { PRIME_MMGBSA }            from '../modules/schrodinger/prime_mmgbsa'

process MERGE_DOCK_CHUNKS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(chunk_sdfs)

    output:
    tuple val(meta), path("${meta.id}_gnina_docked_all.sdf"), emit: merged_sdf

    script:
    """
    cat ${chunk_sdfs} > ${meta.id}_gnina_docked_all.sdf
    """

    stub:
    """
    touch ${meta.id}_gnina_docked_all.sdf
    """
}

process MERGE_GLIDE_CHUNKS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(pv_files), path(csv_files), path(sdf_files)

    output:
    tuple val(meta), path("${meta.id}_dock_pv.maegz"), emit: pv_maegz
    tuple val(meta), path("${meta.id}_dock.csv"), emit: scores_csv
    tuple val(meta), path("${meta.id}_dock.sdf"), emit: poses_sdf

    script:
    """
    # 1. Merge PV files (Schrodinger structcat)
    schrodinger structcat ${pv_files} -omae ${meta.id}_dock_pv.maegz

    # 2. Merge SDF pose files
    cat ${sdf_files} > ${meta.id}_dock.sdf

    # 3. Merge CSV scores with header preservation
    python3 -c "
import sys, csv, os
header = None
rows = []
for f in sys.argv[1:]:
    if os.path.exists(f) and os.path.getsize(f) > 0:
        with open(f, 'r', errors='ignore') as fp:
            r = csv.reader(fp)
            for i, row in enumerate(r):
                if i == 0:
                    if header is None:
                        header = row
                else:
                    rows.append(row)
with open('${meta.id}_dock.csv', 'w', newline='') as out:
    w = csv.writer(out)
    if header:
        w.writerow(header)
    w.writerows(rows)
" ${csv_files}
    """

    stub:
    """
    touch ${meta.id}_dock_pv.maegz
    touch ${meta.id}_dock.csv
    touch ${meta.id}_dock.sdf
    """
}

process MERGE_RESC_CHUNKS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(chunk_sdfs)

    output:
    tuple val(meta), path("${meta.id}_rescored.sdf"), emit: rescored_sdf

    script:
    """
    cat ${chunk_sdfs} > ${meta.id}_rescored.sdf
    """

    stub:
    """
    touch ${meta.id}_rescored.sdf
    """
}

process MERGE_MMGBSA_CHUNKS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(csv_files)

    output:
    tuple val(meta), path("${meta.id}_mmgbsa-out.csv"), emit: scores_csv

    script:
    """
    python3 -c "
import sys, csv, os
header = None
rows = []
for f in sorted(sys.argv[1:]):
    if os.path.exists(f) and os.path.getsize(f) > 0:
        with open(f, 'r', errors='ignore') as fp:
            r = csv.reader(fp)
            for i, row in enumerate(r):
                if i == 0:
                    if header is None:
                        header = row
                else:
                    rows.append(row)
with open('${meta.id}_mmgbsa-out.csv', 'w', newline='') as out:
    w = csv.writer(out)
    if header:
        w.writerow(header)
    w.writerows(rows)
" ${csv_files}
    """

    stub:
    """
    touch ${meta.id}_mmgbsa-out.csv
    """
}

workflow DOCK_LIBRARY {
    take:
    ch_receptor_pdb     // channel: [ val(meta), path(receptor_pdb) ]
    ch_ref_lig_sdf      // channel: [ val(meta), path(ref_lig_sdf) ]
    ch_chunks           // channel: [ val(meta), path(chunk_sdf), val(chunk_idx) ]
    ch_receptor_mae     // channel: [ val(meta), path(receptor_mae) ]
    ch_ref_lig_mae      // channel: [ val(meta), path(ref_lig_mae) ]
    ch_prepared_mae     // channel: [ val(meta), path(prepared_chunk_mae) ]
    dock_engine         // val: "both", "gnina", "glide"
    exhaustiveness      // val: 16
    cnn_scoring         // val: "rescore"
    num_modes           // val: 3
    autobox_add         // val: 4.0
    top_k_frames        // val: 3
    percentile_cutoff   // val: 0.50
    membrane_enabled    // val: true/false

    main:
    ch_summary_files = Channel.empty()
    ch_primary_poses = Channel.empty()

    // 1. GLIDE DOCKING & GNINA RESCORE & SCREENING FUNNEL -> PRIME MM-GBSA
    if (dock_engine in ["glide", "both"]) {
        def inner_box = params.inner_box ?: "10, 10, 10"
        def outer_box = params.outer_box ?: "26, 26, 26"
        def ff_dock   = params.glide_forcefield ?: "OPLS4"
        def prec      = params.glide_precision ?: "SP"
        def poses_n   = params.poses_per_lig ?: 1
        def cnn_resc  = params.cnn_model ?: "crossdock_default2018"
        def prime_job = params.prime_job_type ?: "ENERGY"

        ch_grid_in = ch_receptor_mae.join(ch_ref_lig_mae, by: 0)
        GLIDE_GRIDGEN(ch_grid_in, inner_box, outer_box, ff_dock)

        ch_glide_in = GLIDE_GRIDGEN.out.grid_zip.combine(ch_prepared_mae)
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
        GLIDE_DOCK(ch_glide_in, prec, poses_n, ff_dock)

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
            .groupTuple(by: 0)

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
            .groupTuple(by: 0)

        MERGE_RESC_CHUNKS(ch_to_merge_resc)

        // Collect all consolidated Glide CSVs and GNINA Rescored SDFs across all frames
        ch_all_glide_csvs = MERGE_GLIDE_CHUNKS.out.scores_csv.map { meta, csv -> csv }.collect()
        ch_all_resc_sdfs  = MERGE_RESC_CHUNKS.out.rescored_sdf.map { meta, sdf -> sdf }.collect()

        def bench_ligands = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""

        // Build consensus triage whitelist (Top 50% ligands, Top 3 frames per ligand)
        BUILD_FUNNEL_WHITELIST(
            params.project_name,
            ch_all_glide_csvs,
            ch_all_resc_sdfs,
            top_k_frames,
            percentile_cutoff,
            bench_ligands
        )

        // Filter consolidated PV files per frame and split into parallel chunks
        ch_filter_in = MERGE_GLIDE_CHUNKS.out.pv_maegz
            .combine(BUILD_FUNNEL_WHITELIST.out.whitelist_csv)
            .map { meta, pv, wl -> [ meta, pv ] }

        def batch_pv = params.mmgbsa_batch_size ?: 25
        FILTER_PV(
            ch_filter_in,
            BUILD_FUNNEL_WHITELIST.out.whitelist_csv,
            batch_pv
        )

        // 2D PARALLELISM: N frames * K PV chunks = N * K parallel Prime MM-GBSA Slurm jobs!
        ch_mmgbsa_in = FILTER_PV.out.filtered_pv_chunks.flatMap { meta, files ->
            def fileList = files instanceof List ? files : [files]
            fileList.sort { it.name }.withIndex().collect { f, idx ->
                def chunk_meta = [
                    id: "${meta.id}_mmc${idx}",
                    frame_id: meta.id,
                    frame: meta.frame,
                    chunk_idx: idx,
                    project: meta.project
                ]
                [ chunk_meta, f ]
            }
        }

        PRIME_MMGBSA(ch_mmgbsa_in, "ENERGY", membrane_enabled)

        // Consolidate MM-GBSA chunk CSVs per frame
        ch_mmgbsa_to_merge = PRIME_MMGBSA.out.scores_csv
            .map { chunk_meta, csv ->
                def frame_meta = [
                    id: chunk_meta.frame_id,
                    frame: chunk_meta.frame,
                    project: chunk_meta.project
                ]
                [ frame_meta, csv ]
            }
            .groupTuple(by: 0)

        MERGE_MMGBSA_CHUNKS(ch_mmgbsa_to_merge)

        ch_summary_files = ch_summary_files
            .mix(MERGE_GLIDE_CHUNKS.out.scores_csv.map { meta, csv -> csv })
            .mix(MERGE_GLIDE_CHUNKS.out.poses_sdf.map { meta, sdf -> sdf })
            .mix(MERGE_RESC_CHUNKS.out.rescored_sdf.map { meta, sdf -> sdf })
            .mix(BUILD_FUNNEL_WHITELIST.out.whitelist_csv)
            .mix(MERGE_MMGBSA_CHUNKS.out.scores_csv.map { meta, csv -> csv })

        ch_primary_poses = MERGE_GLIDE_CHUNKS.out.poses_sdf
    }

    // 2. GNINA DOCKING
    if (dock_engine in ["gnina", "both"]) {
        ch_rec_comb = ch_receptor_pdb.join(ch_ref_lig_sdf, by: 0)
        ch_dock_in = ch_chunks.combine(ch_rec_comb)
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
            .groupTuple(by: 0)

        MERGE_DOCK_CHUNKS(ch_to_merge)

        ch_summary_files = ch_summary_files
            .mix(MERGE_DOCK_CHUNKS.out.merged_sdf.map { meta, sdf -> sdf })

        if (dock_engine == "gnina") {
            ch_primary_poses = MERGE_DOCK_CHUNKS.out.merged_sdf
        }
    }

    // Collect all summary files across all frames into a single flat list
    ch_all_summary_files = ch_summary_files.collect()

    emit:
    docked_poses   = ch_primary_poses
    summary_files  = ch_all_summary_files
}
