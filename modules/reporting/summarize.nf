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
    path(summary_files, stageAs: 'in_files/*')
    val(native_ligand)
    val(benchmark_ligands)

    output:
    path("ranking.csv"), emit: ranking_csv
    path("cross_frame_matrix.csv"), emit: matrix_csv
    path("summary.json"), emit: summary_json

    script:
    def nat_arg = native_ligand ? "--native_ligand ${native_ligand}" : ""
    def bench_arg = benchmark_ligands ? "--benchmark_ligands ${benchmark_ligands}" : ""
    """
    cp -L -t . in_files/*

    PLAN_ARG=""
    if [ -f rescore_plan.json ]; then
        PLAN_ARG="--plan rescore_plan.json"
    elif [ -f dock_plan.json ]; then
        PLAN_ARG="--plan dock_plan.json"
    fi

    python3 ${projectDir}/bin/summarize_results.py \\
        --docking_dir . \\
        \${PLAN_ARG} \\
        ${nat_arg} \\
        ${bench_arg} \\
        --outdir .
    """

    stub:
    """
    touch ranking.csv
    touch cross_frame_matrix.csv
    touch summary.json
    """
}
