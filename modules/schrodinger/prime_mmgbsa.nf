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
    def mem_flag = (membrane_enabled.toString() == 'true') ? "-membrane" : ""
    """
    schrodinger prime_mmgbsa \
        -WAIT \
        -LOCAL \
        -job_type ${job_type} \
        -csv_output yes \
        -out_type PV \
        -jobname ${meta.id}_mmgbsa \
        -OVERWRITE \
        ${mem_flag} \
        ${pv_maegz}

    # Normalize output names if Schrödinger generated without -out infix
    if [ -f "${meta.id}_mmgbsa.csv" ] && [ ! -f "${meta.id}_mmgbsa-out.csv" ]; then
        mv "${meta.id}_mmgbsa.csv" "${meta.id}_mmgbsa-out.csv"
    fi
    if [ -f "${meta.id}_mmgbsa.maegz" ] && [ ! -f "${meta.id}_mmgbsa-out.maegz" ]; then
        mv "${meta.id}_mmgbsa.maegz" "${meta.id}_mmgbsa-out.maegz"
    fi

    if [ ! -f "${meta.id}_mmgbsa-out.csv" ] || [ ! -f "${meta.id}_mmgbsa-out.maegz" ]; then
        echo "FATAL: Prime MM-GBSA failed to generate expected outputs."
        echo "Directory listing:"
        ls -la
        if [ -f "${meta.id}_mmgbsa.log" ]; then
            echo "=== ${meta.id}_mmgbsa.log ==="
            cat "${meta.id}_mmgbsa.log"
        fi
        exit 1
    fi    """
    stub:
    """
    touch ${meta.id}_mmgbsa-out.csv
    touch ${meta.id}_mmgbsa-out.maegz
    """
}