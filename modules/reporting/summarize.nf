/*
 * modules/reporting/summarize.nf
 * Aggregates all raw docking, MM-GBSA, and GNINA metrics into consensus ranking CSV tables
 */

process SUMMARIZE_DOCKING {
    tag "${project_name}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val(project_name)
    path(summary_files)
    val(native_ligand)
    val(benchmark_ligands)

    output:
    path("ranking.csv"), emit: ranking_csv
    path("cross_frame_matrix.csv"), emit: matrix_csv
    path("summary.json"), emit: summary_json

    script:
    """
    summarize_results.py \\
        --docking_dir . \\
        --outdir .
    """

    stub:
    """
    touch ranking.csv
    touch cross_frame_matrix.csv
    touch summary.json
    """
}
