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
    def lig_title = meta.ligand_title ?: meta.id
    """
    IN_SDF=\$(ls in_ligand/* | head -n 1)
    if [ -z "\$IN_SDF" ] || [ ! -s "\$IN_SDF" ]; then
        echo "FATAL: Ligand file is missing or empty in CONVERT_READY_LIGAND!"
        exit 1
    fi

    schrodinger run python3 -c "
from schrodinger import structure
st = structure.StructureReader.read('\$IN_SDF')
title = '${lig_title}'.strip()
st.title = title
st.property['s_m_title'] = title
st.write('${meta.id}.maegz')
st.write('${meta.id}.sdf')
"
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
            def chunk_meta = [
                id: "${meta.id}_c${idx}",
                name: meta.name ?: meta.id,
                chunk_idx: idx,
                project: meta.project ?: meta.id,
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
                def chunk_meta = [
                    id: "${meta.id}_ready_${idx}",
                    name: meta.name ?: meta.id,
                    chunk_idx: "ready_${idx}",
                    project: meta.project ?: meta.id,
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
