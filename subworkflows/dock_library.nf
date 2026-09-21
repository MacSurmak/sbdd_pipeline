/*
 * subworkflows/dock_library.nf
 * Distributed batch docking across Slurm GPU worker nodes
 */

include { GNINA_DOCK } from '../modules/gnina/dock'
include { GLIDE_DOCK } from '../modules/schrodinger/glide_dock'

process MERGE_DOCK_CHUNKS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(chunk_sdfs)

    output:
    tuple val(meta), path("${meta.id}_docked_all.sdf"), emit: merged_sdf

    script:
    """
    cat chunk_sdfs* > ${meta.id}_docked_all.sdf 2>/dev/null || cat *.sdf > ${meta.id}_docked_all.sdf
    cat ${chunk_sdfs} > ${meta.id}_docked_all.sdf
    """

    stub:
    """
    touch ${meta.id}_docked_all.sdf
    """
}

workflow DOCK_LIBRARY {
    take:
    ch_receptor_pdb     // channel: [ val(meta), path(receptor_pdb) ]
    ch_ref_lig_sdf      // channel: [ val(meta), path(ref_lig_sdf) ]
    ch_chunks           // channel: [ val(meta), path(chunk_sdf), val(chunk_idx) ]
    dock_engine         // val: "gnina" or "glide"
    exhaustiveness      // val: 16
    cnn_scoring         // val: "rescore"
    num_modes           // val: 3
    autobox_add         // val: 4.0

    main:
    if (dock_engine == "gnina") {
        // Combine inputs for GNINA_DOCK: [ meta, receptor_pdb, ref_lig_sdf, chunk_sdf, chunk_idx ]
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

        // Group all chunk outputs by meta to merge
        ch_to_merge = GNINA_DOCK.out.docked_sdf
            .groupTuple(by: 0)

        MERGE_DOCK_CHUNKS(ch_to_merge)
        ch_final_poses = MERGE_DOCK_CHUNKS.out.merged_sdf
        ch_logs = GNINA_DOCK.out.log
    } else {
        ch_final_poses = Channel.empty()
        ch_logs = Channel.empty()
    }

    emit:
    docked_poses = ch_final_poses
    logs         = ch_logs
}

