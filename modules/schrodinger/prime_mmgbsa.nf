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
    schrodinger prime_mmgbsa \
        -job_type ${job_type} \
        -csv_output yes \
        -out_type PV \
        -jobname ${meta.id}_mmgbsa \
        -LOCAL \
        -NJOBS ${task.cpus} \
        -OVERWRITE \
        ${mem_flag} \
        ${pv_maegz}

    # Defensive check: ensure output files exist even if a chunk had zero scorable ligands
    if [ ! -f "${meta.id}_mmgbsa-out.csv" ]; then
        echo "title,r_psp_MMGBSA_dG_Bind" > "${meta.id}_mmgbsa-out.csv"
    fi
    if [ ! -f "${meta.id}_mmgbsa-out.maegz" ]; then
        cp "${pv_maegz}" "${meta.id}_mmgbsa-out.maegz"
    fi    """

    stub:
    """
    touch ${meta.id}_mmgbsa-out.csv
    touch ${meta.id}_mmgbsa-out.maegz
    """
}
