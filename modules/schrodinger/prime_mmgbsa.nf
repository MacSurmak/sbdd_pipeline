/*
 * modules/schrodinger/prime_mmgbsa.nf
 * Prime MM-GBSA Free Energy Calculation
 */

process PRIME_MMGBSA {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(pv_maegz)
    val(job_type)       // ENERGY or REAL_MIN
    val(membrane_enabled) // boolean

    output:
    tuple val(meta), path("${meta.id}_mmgbsa-out.csv"), emit: scores_csv
    tuple val(meta), path("${meta.id}_mmgbsa-out.maegz"), emit: out_maegz

    script:
    def mem_flag = (membrane_enabled.toString() == 'true') ? "-membrane" : ""
    """
    schrodinger prime_mmgbsa \\
        -job_type ${job_type} \\
        -csv_output yes \\
        -out_type PV \\
        -jobname ${meta.id}_mmgbsa \\
        -HOST localhost:${task.cpus} \\
        -NJOBS ${task.cpus} \\
        -OVERWRITE \\
        -WAIT \\
        ${mem_flag} \\
        ${pv_maegz}
    """

    stub:
    """
    touch ${meta.id}_mmgbsa-out.csv
    touch ${meta.id}_mmgbsa-out.maegz
    """
}
