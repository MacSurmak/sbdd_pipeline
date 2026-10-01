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
    tuple val(meta), path(sdf, stageAs: 'ligand.sdf'), val(idx)

    output:
    tuple val(meta), path("${meta.id}.maegz"), emit: prepared_maegz
    tuple val(meta), path(sdf), emit: prepared_sdf

    script:
    """
    if [ ! -s "ligand.sdf" ]; then
        echo "FATAL: Ligand file ${sdf} is empty! The ligand was likely not found in the receptor during extraction."
        exit 1
    fi

    schrodinger run python3 -c "
import sys
from schrodinger import structure
try:
    with structure.StructureReader('ligand.sdf') as reader, structure.StructureWriter('${meta.id}.maegz') as writer:
        for st in reader:
            writer.append(st)
except Exception as e:
    print(f'Conversion failed: {e}', file=sys.stderr)
    sys.exit(1)
"
    """

    stub:
    """
    touch ${meta.id}.maegz
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

    // 3. Ready libraries bypass LigPrep entirely
    ch_ready_in = branched_libs.ready.map { meta, f ->
        def chunk_meta = [
            id: "${meta.id}_ready",
            chunk_idx: "native",
            project: meta.project ?: meta.id,
            status: "ready",
            is_native: meta.is_native ?: false
        ]
        [ chunk_meta, f, "native" ]
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
