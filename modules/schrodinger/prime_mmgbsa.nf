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
    # Dynamically bind full host RPC, keyutils and Kerberos library chain into container
    for lib in libtirpc.so.3 libkeyutils.so.1 libgssapi_krb5.so.2 libkrb5.so.3 libk5crypto.so.3 libcom_err.so.2 libkrb5support.so.0; do
        LIB_PATH=\$(ls /lib/x86_64-linux-gnu/\$lib /usr/lib/x86_64-linux-gnu/\$lib /usr/lib/\$lib 2>/dev/null | head -n 1)
        if [ -n "\$LIB_PATH" ]; then
            export APPTAINER_BINDPATH="\${LIB_PATH}:/usr/lib/x86_64-linux-gnu/\$lib,\${LIB_PATH}:/lib/x86_64-linux-gnu/\$lib\${APPTAINER_BINDPATH:+,\$APPTAINER_BINDPATH}"
            export SINGULARITY_BINDPATH="\${LIB_PATH}:/usr/lib/x86_64-linux-gnu/\$lib,\${LIB_PATH}:/lib/x86_64-linux-gnu/\$lib\${SINGULARITY_BINDPATH:+,\$SINGULARITY_BINDPATH}"
        fi
    done

    export SCHRODINGER_JOBDB="/tmp/schrodinger_jobdb_\${SLURM_JOB_ID:-\$\$}"
    mkdir -p "\$SCHRODINGER_JOBDB"

    schrodinger prime_mmgbsa \
        -job_type ${job_type} \
        -csv_output yes \
        -out_type PV \
        -jobname ${meta.id}_mmgbsa \
        -HOST localhost \
        -NJOBS 1 \
        -OVERWRITE \
        -WAIT \
        ${mem_flag} \
        ${pv_maegz}

    # Normalize CSV name: Prime MM-GBSA writes <jobname>.csv, not <jobname>-out.csv
    if [ -f "${meta.id}_mmgbsa.csv" ]; then
        cp "${meta.id}_mmgbsa.csv" "${meta.id}_mmgbsa-out.csv"
    fi
    if [ -f "${meta.id}_mmgbsa.maegz" ] && [ ! -f "${meta.id}_mmgbsa-out.maegz" ]; then
        mv "${meta.id}_mmgbsa.maegz" "${meta.id}_mmgbsa-out.maegz"
    fi
    if [ ! -f "${meta.id}_mmgbsa-out.csv" ] || [ ! -f "${meta.id}_mmgbsa-out.maegz" ]; then
        echo "FATAL: Prime MM-GBSA failed to generate expected outputs."
        echo "Directory listing:"
        ls -la
        if [ -f "${meta.id}_mmgbsa.Prime.log" ]; then
            echo "=== ${meta.id}_mmgbsa.Prime.log ==="
            cat "${meta.id}_mmgbsa.Prime.log"
        elif [ -f "${meta.id}_mmgbsa.log" ]; then
            echo "=== ${meta.id}_mmgbsa.log ==="
            cat "${meta.id}_mmgbsa.log"
        fi
        exit 1
    fi
    """

    stub:
    """
    touch ${meta.id}_mmgbsa-out.csv
    touch ${meta.id}_mmgbsa-out.maegz
    """
}