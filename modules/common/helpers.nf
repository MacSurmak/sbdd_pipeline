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

process BATCH_DOCK_LIGANDS {
    tag "${project_name}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val(project_name)
    path(dock_plan_json)
    val(batch_size)

    output:
    path("chunks_manifest.json"), emit: manifest
    path("chunk_*.sdf"), emit: sdf_chunks
    path("chunk_*.maegz"), emit: mae_chunks

    script:
    """
    python3 ${projectDir}/bin/batch_dock_ligands.py \\
        --dock_plan ${dock_plan_json} \\
        --batch_size ${batch_size} \\
        --out_manifest chunks_manifest.json
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
