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
}

process BATCH_UNIDOCK_LIGANDS {
    tag "${project_name}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val(project_name)
    path(dock_plan_json)
    val(unidock_batch_size)

    output:
    path("unidock_chunks_manifest.json"), emit: manifest
    path("unidock_chunk_*.sdf"), emit: sdf_chunks

    script:
    """
    python3 ${projectDir}/bin/batch_dock_ligands.py \\
        --dock_plan ${dock_plan_json} \\
        --glide_batch_size 0 \\
        --unidock_batch_size ${unidock_batch_size} \\
        --out_manifest unidock_chunks_manifest.json
    """
}

process BATCH_GLIDE_LIGANDS {
    tag "${project_name}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val(project_name)
    path(dock_plan_json)
    val(glide_batch_size)

    output:
    path("glide_chunks_manifest.json"), emit: manifest
    path("glide_chunk_*.maegz"), emit: mae_chunks

    script:
    """
    python3 ${projectDir}/bin/batch_dock_ligands.py \\
        --dock_plan ${dock_plan_json} \\
        --glide_batch_size ${glide_batch_size} \\
        --unidock_batch_size 0 \\
        --out_manifest glide_chunks_manifest.json
    """
}

process MERGE_CSV {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(csv_files, stageAs: 'csvs_?/*')

    output:
    tuple val(meta), path("merged.csv"), emit: merged_csv

    script:
    """
    python3 ${projectDir}/bin/merge_csv.py --inputs ${csv_files} --output merged.csv
    """
}
