/*
 * subworkflows/prepare_ligands.nf
 * Splits SDF library into parallel chunks and/or runs Schrodinger LigPrep
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
    dock_engine   // val: "both", "gnina", "glide"
    batch_size    // val: 15
    target_ph     // val: 7.4
    ph_tol        // val: 2.0
    ff_version    // val: 16
    max_states    // val: 32

    main:
    ch_prepared_mae = Channel.empty()
    ch_chunks       = Channel.empty()

    if (dock_engine in ["glide", "both"]) {
        LIGPREP(ch_library, target_ph, ph_tol, ff_version, max_states)
        ch_prepared_mae = LIGPREP.out.prepared_maegz
    }

    if (dock_engine in ["gnina", "both"]) {
        SPLIT_SDF(ch_library, batch_size)
        ch_chunks = SPLIT_SDF.out.chunks.flatMap { meta, files ->
            def fileList = files instanceof List ? files : [files]
            fileList.withIndex().collect { f, idx -> [ meta, f, idx ] }
        }
    }

    emit:
    prepared_ligands = ch_prepared_mae
    chunks           = ch_chunks
}
