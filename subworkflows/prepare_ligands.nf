/*
 * subworkflows/prepare_ligands.nf
 * Splits SDF library into parallel chunks and/or runs Schrodinger LigPrep
 */

include { LIGPREP } from '../modules/schrodinger/ligprep'
include { SPLIT_SDF } from '../modules/common/helpers'
include { CONVERT_READY_LIGAND } from '../modules/schrodinger/helpers'

workflow PREPARE_LIGANDS {
    take:
    ch_library    // channel: [ val(meta), path(library_file) ]
    batch_size    // val: 15
    target_ph     // val: 7.4
    ph_tol        // val: 2.0
    ff_version    // val: 16
    max_states    // val: 32

    main:
    // 1. Branch libraries based on status
    ch_library.branch {
        ready: it[0].status == 'ready'
        prep:  true
    }.set { branched_libs }

    // 2. Unprepared ligands/libraries undergo Splitting and LigPrep
    SPLIT_SDF(branched_libs.prep, batch_size)
    ch_prep_chunks = SPLIT_SDF.out.chunks.flatMap { meta, files ->
        def fileList = files instanceof List ? files : (files ? [files] : [])
        fileList.withIndex().collect { f, idx ->
            def chunk_meta = meta + [
                id: "${meta.id}_c${idx}",
                chunk_idx: idx,
                status: 'ready'
            ]
            [ chunk_meta, f ]
        }
    }

    LIGPREP(ch_prep_chunks, target_ph, ph_tol, ff_version, max_states)

    // 3. Ready libraries bypass LigPrep
    ch_ready_in = branched_libs.ready
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .flatMap { list ->
            list.withIndex().collect { item, idx ->
                def meta = item[0]
                def f = item[1]
                def chunk_meta = meta + [
                    id: "${meta.id}_ready_${idx}",
                    chunk_idx: "ready_${idx}",
                    status: "ready"
                ]
                [ chunk_meta, f, "ready_${idx}" ]
            }
        }
    CONVERT_READY_LIGAND(ch_ready_in)

    // 4. Combine outputs
    ch_prepared_mae = LIGPREP.out.prepared_maegz.mix(CONVERT_READY_LIGAND.out.prepared_maegz)
    ch_prepared_sdf = LIGPREP.out.prepared_sdf.mix(CONVERT_READY_LIGAND.out.prepared_sdf)

    emit:
    prepared_ligands = ch_prepared_mae
    prepared_sdfs    = ch_prepared_sdf
}
