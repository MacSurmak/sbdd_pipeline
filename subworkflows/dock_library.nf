/*
 * subworkflows/dock_library.nf
 * Multi-Engine Docking & Rescoring Pipeline (Glide, Prime MM-GBSA, GNINA, CNN Rescore)
 */

include { GNINA_DOCK }    from '../modules/gnina/dock'
include { GNINA_RESCORE } from '../modules/gnina/rescore'
include { GLIDE_GRIDGEN } from '../modules/schrodinger/gridgen'
include { GLIDE_DOCK }    from '../modules/schrodinger/glide_dock'
include { PRIME_MMGBSA }  from '../modules/schrodinger/prime_mmgbsa'

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
    cat chunk_sdfs* > ${meta.id}_gnina_docked_all.sdf 2>/dev/null || cat *.sdf > ${meta.id}_gnina_docked_all.sdf
    cat ${chunk_sdfs} > ${meta.id}_gnina_docked_all.sdf
    """

    stub:
    """
    touch ${meta.id}_gnina_docked_all.sdf
    """
}

workflow DOCK_LIBRARY {
    take:
    ch_receptor_pdb     // channel: [ val(meta), path(receptor_pdb) ]
    ch_ref_lig_sdf      // channel: [ val(meta), path(ref_lig_sdf) ]
    ch_chunks           // channel: [ val(meta), path(chunk_sdf), val(chunk_idx) ]
    ch_receptor_mae     // channel: [ val(meta), path(receptor_mae) ]
    ch_ref_lig_mae      // channel: [ val(meta), path(ref_lig_mae) ]
    ch_prepared_mae     // channel: [ val(meta), path(prepared_mae) ]
    dock_engine         // val: "both", "gnina", "glide"
    exhaustiveness      // val: 16
    cnn_scoring         // val: "rescore"
    num_modes           // val: 3
    autobox_add         // val: 4.0

    main:
    ch_summary_files = Channel.empty()
    ch_primary_poses = Channel.empty()

    // 1. GLIDE DOCKING & PRIME MM-GBSA & GNINA RESCORE
    if (dock_engine in ["glide", "both"]) {
        ch_grid_in = ch_receptor_mae.join(ch_ref_lig_mae, by: 0)
        GLIDE_GRIDGEN(ch_grid_in, "10, 10, 10", "26, 26, 26", "OPLS4")

        ch_glide_in = GLIDE_GRIDGEN.out.grid_zip.join(ch_prepared_mae, by: 0)
        GLIDE_DOCK(ch_glide_in, "SP", 1, "OPLS4")

        PRIME_MMGBSA(GLIDE_DOCK.out.pv_maegz, "ENERGY")

        ch_rescore_in = ch_receptor_pdb.join(GLIDE_DOCK.out.poses_sdf, by: 0)
        GNINA_RESCORE(ch_rescore_in, "crossdock_default2018")

        ch_summary_files = ch_summary_files
            .mix(GLIDE_DOCK.out.scores_csv.map { meta, csv -> [ meta, [csv] ] })
            .mix(PRIME_MMGBSA.out.scores_csv.map { meta, csv -> [ meta, [csv] ] })
            .mix(GNINA_RESCORE.out.rescored_sdf.map { meta, sdf -> [ meta, [sdf] ] })

        ch_primary_poses = GLIDE_DOCK.out.poses_sdf
    }

    // 2. GNINA DOCKING
    if (dock_engine in ["gnina", "both"]) {
        ch_dock_in = ch_chunks.combine(ch_receptor_pdb, by: 0).combine(ch_ref_lig_sdf, by: 0)
            .map { meta, chunk_sdf, chunk_idx, rec_pdb, ref_sdf ->
                [ meta, rec_pdb, ref_sdf, chunk_sdf, chunk_idx ]
            }

        GNINA_DOCK(
            ch_dock_in,
            exhaustiveness,
            cnn_scoring,
            num_modes,
            autobox_add
        )

        ch_to_merge = GNINA_DOCK.out.docked_sdf.groupTuple(by: 0)
        MERGE_DOCK_CHUNKS(ch_to_merge)

        ch_summary_files = ch_summary_files
            .mix(MERGE_DOCK_CHUNKS.out.merged_sdf.map { meta, sdf -> [ meta, [sdf] ] })

        if (dock_engine == "gnina") {
            ch_primary_poses = MERGE_DOCK_CHUNKS.out.merged_sdf
        }
    }

    // Group all summary files (SDFs, CSVs) by meta for master ranking
    ch_bundled_for_summary = ch_summary_files
        .groupTuple(by: 0)
        .map { meta, file_lists ->
            def all_files = file_lists.flatten()
            [ meta, all_files ]
        }

    emit:
    docked_poses   = ch_primary_poses
    summary_bundle = ch_bundled_for_summary
}
