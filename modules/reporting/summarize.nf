/*
 * modules/reporting/summarize.nf
 * Aggregates all raw docking, MM-GBSA, and MD metrics into ranking tables
 */

process SUMMARIZE_DOCKING {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(docked_poses)

    output:
    tuple val(meta), path("ranking.csv"), emit: ranking_csv
    tuple val(meta), path("summary.json"), emit: summary_json

    script:
    """
    python3 ${projectDir}/bin/summarize_results.py \\
        --docking_dir . \\
        --outdir .
    """

    stub:
    """
    touch ranking.csv
    touch summary.json
    """
}
