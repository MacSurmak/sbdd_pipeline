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
