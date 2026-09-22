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

process MERGE_LIGPREP {
    tag "${project_name}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val(project_name)
    path(chunk_maes)

    output:
    tuple val([id: project_name, project: project_name]), path("${project_name}_prepared.maegz"), emit: prepared_maegz

    script:
    """
    schrodinger structcat ${chunk_maes} -omae ${project_name}_prepared.maegz
    """

    stub:
    """
    touch ${project_name}_prepared.maegz
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

    // 1. Split SDF library into chunks for distributed computing across Slurm
    SPLIT_SDF(ch_library, batch_size)
    ch_chunks = SPLIT_SDF.out.chunks.flatMap { meta, files ->
        def fileList = files instanceof List ? files : [files]
        fileList.withIndex().collect { f, idx ->
            def chunk_meta = [
                id: "${meta.id}_c${idx}",
                chunk_idx: idx,
                project: meta.project ?: meta.id
            ]
            [ chunk_meta, f, idx ]
        }
    }

    // 2. Parallel LIGPREP across Slurm worker nodes
    if (dock_engine in ["glide", "both"]) {
        ch_ligprep_in = ch_chunks.map { chunk_meta, f, idx -> [ chunk_meta, f ] }
        LIGPREP(ch_ligprep_in, target_ph, ph_tol, ff_version, max_states)

        ch_all_maegz = LIGPREP.out.prepared_maegz
            .map { chunk_meta, maegz -> maegz }
            .collect()

        MERGE_LIGPREP(params.project_name, ch_all_maegz)
        ch_prepared_mae = LIGPREP.out.prepared_maegz
    }

    emit:
    prepared_ligands = ch_prepared_mae
    chunks           = ch_chunks
}
