/*
 * modules/schrodinger/prime_mmgbsa.nf
 * Prime MM-GBSA Free Energy Calculation
 */
process PRIME_MMGBSA {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(pv_maegz)
    val(job_type)       // ENERGY or REAL_MIN
    val(membrane_enabled) // boolean

    output:
    tuple val(meta), path("${meta.id}_mmgbsa-out.csv"), emit: scores_csv
    tuple val(meta), path("${meta.id}_mmgbsa-out.maegz"), emit: out_maegz

    script:
    def mem_arg = (membrane_enabled.toString() == 'true') ? "--membrane" : ""
    """
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py prime_mmgbsa \\
        --pv_maegz ${pv_maegz} \\
        --job_type ${job_type} \\
        --job_id ${meta.id} \\
        ${mem_arg}
    """
}