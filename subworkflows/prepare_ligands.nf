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
    ch_ref_native // channel: [ val(meta), path(ref_native_sdf) ]
    dock_engine   // val: "both", "gnina", "glide"
    batch_size    // val: 15
    target_ph     // val: 7.4
    ph_tol        // val: 2.0
    ff_version    // val: 16
    max_states    // val: 32

    main:
    // 1. Split raw SDF library into chunks for distributed computing across Slurm
    SPLIT_SDF(ch_library, batch_size)
    ch_raw_chunks = SPLIT_SDF.out.chunks.flatMap { meta, files ->
        def fileList = files instanceof List ? files : [files]
        fileList.withIndex().collect { f, idx ->
            def chunk_meta = [
                id: "${meta.id}_c${idx}",
                chunk_idx: idx,
                project: meta.project ?: meta.id,
                is_native: false
            ]
            [ chunk_meta, f, idx ]
        }
    }

    // 2. Add native reference ligand as a distinct chunk so it is screened and rescored in all frames
    ch_native_chunk = ch_ref_native.first().map { meta, native_sdf ->
        def nat_meta = [
            id: "${meta.project ?: params.project_name}_native",
            chunk_idx: 9999,
            project: meta.project ?: params.project_name,
            is_native: true
        ]
        [ nat_meta, native_sdf, 9999 ]
    }

    ch_all_chunks = ch_raw_chunks.mix(ch_native_chunk)

    // 3. Prepare 3D protonated / Epik states at physiological pH across Slurm
    ch_ligprep_in = ch_all_chunks.map { chunk_meta, f, idx -> [ chunk_meta, f ] }
    LIGPREP(ch_ligprep_in, target_ph, ph_tol, ff_version, max_states)

    ch_prepared_mae = LIGPREP.out.prepared_maegz
    ch_prepared_sdf = LIGPREP.out.prepared_sdf.map { chunk_meta, sdf ->
        [ chunk_meta, sdf, chunk_meta.chunk_idx ]
    }

    ch_all_maegz = LIGPREP.out.prepared_maegz
        .map { chunk_meta, maegz -> maegz }
        .collect()

    MERGE_LIGPREP(params.project_name, ch_all_maegz)

    emit:
    prepared_ligands = ch_prepared_mae
    chunks           = ch_prepared_sdf
    merged_library   = MERGE_LIGPREP.out.prepared_maegz
}
