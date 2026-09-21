/*
 * subworkflows/prepare_ligands.nf
 * Splits SDF library into parallel chunks or runs Schrodinger LigPrep
 */

include { LIGPREP } from '../modules/schrodinger/ligprep'

process SPLIT_SDF {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(library_sdf)
    val(batch_size)

    output:
    tuple val(meta), path("chunk_*.sdf"), emit: chunks

    script:
    """
    python3 ${projectDir}/bin/split_sdf.py \\
        --input ${library_sdf} \\
        --batch_size ${batch_size} \\
        --outdir .
    """

    stub:
    """
    touch chunk_000.sdf
    touch chunk_001.sdf
    """
}

workflow PREPARE_LIGANDS {
    take:
    ch_library    // channel: [ val(meta), path(library_file) ]
    dock_engine   // val: "gnina" or "glide"
    batch_size    // val: 15
    target_ph     // val: 7.4
    ph_tol        // val: 2.0
    ff_version    // val: 16
    max_states    // val: 32

    main:
    if (dock_engine == "glide") {
        LIGPREP(ch_library, target_ph, ph_tol, ff_version, max_states)
        ch_prepared = LIGPREP.out.prepared_maegz
        ch_chunks = Channel.empty()
    } else {
        SPLIT_SDF(ch_library, batch_size)
        // Transpose chunks to emit [ meta, chunk_path, chunk_idx ]
        ch_chunks = SPLIT_SDF.out.chunks.flatMap { meta, files ->
            def fileList = files instanceof List ? files : [files]
            fileList.withIndex().collect { f, idx -> [ meta, f, idx ] }
        }
        ch_prepared = ch_library
    }

    emit:
    prepared_ligands = ch_prepared
    chunks           = ch_chunks
}

