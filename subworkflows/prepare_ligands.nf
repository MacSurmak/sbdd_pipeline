/*
 * subworkflows/prepare_ligands.nf
 * Splits SDF library into parallel chunks and/or runs Schrodinger LigPrep
 */

include { LIGPREP   } from '../modules/schrodinger/ligprep'
include { SPLIT_SDF } from '../modules/common/helpers'

workflow PREPARE_LIGANDS {
    take:
    ch_library    // channel: [ val(meta), path(library_file) ]
    batch_size    // val: 15
    target_ph     // val: 7.4
    ph_tol        // val: 2.0
    ff_version    // val: 16
    max_states    // val: 32

    main:
    SPLIT_SDF(ch_library, batch_size)
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

    emit:
    prepared_ligands = LIGPREP.out.prepared_maegz
    prepared_sdfs    = LIGPREP.out.prepared_sdf
}
