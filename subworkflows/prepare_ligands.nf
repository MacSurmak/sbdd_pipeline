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

process CONVERT_READY_LIGAND {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sdf, stageAs: 'in_ligand/*'), val(idx)

    output:
    tuple val(meta), path("${meta.id}.maegz"), emit: prepared_maegz
    tuple val(meta), path("${meta.id}.sdf"), emit: prepared_sdf

    script:
    """
    IN_SDF=\$(ls in_ligand/* | head -n 1)
    if [ ! -s "\$IN_SDF" ]; then
        echo "FATAL: Ligand file \$IN_SDF is empty!"
        exit 1
    fi

    schrodinger structconvert "\$IN_SDF" "${meta.id}.maegz"
    cp "\$IN_SDF" "${meta.id}.sdf"
    """

    stub:
    """
    touch ${meta.id}.maegz
    touch ${meta.id}.sdf
    """
}

process MERGE_LIGPREP {
    tag "${project_name}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val(project_name)
    path(chunk_maes, stageAs: 'chunks_?/*')

    output:
    tuple val([id: project_name, project: project_name]), path("${project_name}_prepared.maegz"), emit: prepared_maegz

    script:
    """
    schrodinger structcat chunks_*/*.mae* -omae ${project_name}_prepared.maegz
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
    // 1. Branch libraries *before* processing based on DB status
    ch_library.branch {
        ready: it[0].status == 'ready'
        prep:  true
    }.set { branched_libs }

    // 2. Unprepared libraries go through splitting and full LigPrep
    SPLIT_SDF(branched_libs.prep, batch_size)
    ch_prep_chunks = SPLIT_SDF.out.chunks.flatMap { meta, files ->
        def fileList = files instanceof List ? files : (files ? [files] : [])
        fileList.withIndex().collect { f, idx ->
            def chunk_meta = [
                id: "${meta.id}_c${idx}",
                chunk_idx: idx,
                project: meta.project ?: meta.id,
                status: meta.status ?: 'unprepared',
                is_native: meta.is_native ?: false
            ]
            [ chunk_meta, f ]
        }
    }

    LIGPREP(ch_prep_chunks, target_ph, ph_tol, ff_version, max_states)

    // 3. Ready libraries bypass LigPrep entirely with unique chunk names
    ch_ready_in = branched_libs.ready
        .toSortedList { a, b -> a[0].id <=> b[0].id }
        .flatMap { list ->
            list.withIndex().collect { item, idx ->
                def meta = item[0]
                def f = item[1]
                def chunk_meta = [
                    id: "${meta.id}_ready_${idx}",
                    chunk_idx: "native_${idx}",
                    project: meta.project ?: meta.id,
                    status: "ready",
                    is_native: meta.is_native ?: false
                ]
                [ chunk_meta, f, "native_${idx}" ]
            }
        }
    CONVERT_READY_LIGAND(ch_ready_in)

    // 4. Mix outputs
    ch_prepared_mae = LIGPREP.out.prepared_maegz.mix(CONVERT_READY_LIGAND.out.prepared_maegz)
    
    ch_prepared_sdf = LIGPREP.out.prepared_sdf
        .mix(CONVERT_READY_LIGAND.out.prepared_sdf)
        .map { chunk_meta, sdf -> [ chunk_meta, sdf, chunk_meta.chunk_idx ] }

    ch_all_maegz = ch_prepared_mae
        .map { chunk_meta, maegz -> maegz }
        .collect()

    MERGE_LIGPREP(params.project_name, ch_all_maegz)

    emit:
    prepared_ligands = ch_prepared_mae
    chunks           = ch_prepared_sdf
    merged_library   = MERGE_LIGPREP.out.prepared_maegz
}
