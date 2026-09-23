/*
 * modules/reporting/build_funnel.nf
 * Builds funnel whitelist CSV using Glide CSVs and GNINA Rescore SDFs
 */

process BUILD_FUNNEL_WHITELIST {
    tag "${project_name}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    val(project_name)
    path(glide_csvs)
    path(gnina_sdfs)
    val(top_k_frames)
    val(percentile_cutoff)
    val(benchmark_ligands)
    val(native_ligand)

    output:
    path("funnel_whitelist.csv"), emit: whitelist_csv

    script:
    def bench_arg = benchmark_ligands ? "--benchmark_ligands ${benchmark_ligands}" : ""
    def nat_arg = native_ligand ? "--native_ligand ${native_ligand}" : ""
    """
    ${params.workdir}/bin/build_funnel_whitelist.py \\
        --glide_csvs ${glide_csvs} \\
        --gnina_sdfs ${gnina_sdfs} \\
        --top_k_frames ${top_k_frames} \\
        --percentile_cutoff ${percentile_cutoff} \\
        ${bench_arg} \\
        ${nat_arg} \\
        --output funnel_whitelist.csv
    """

    stub:
    """
    touch funnel_whitelist.csv
    """
}
