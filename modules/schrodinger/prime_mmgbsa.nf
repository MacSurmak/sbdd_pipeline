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
    # Ensure valid non-empty machine-id for Schrödinger mmlic3 in container
    if [ -s /etc/machine-id ]; then
        export APPTAINER_BINDPATH="/etc/machine-id\${APPTAINER_BINDPATH:+,\$APPTAINER_BINDPATH}"
        export SINGULARITY_BINDPATH="/etc/machine-id\${SINGULARITY_BINDPATH:+,\$SINGULARITY_BINDPATH}"
    else
        FAKE_MID=\$(mktemp /tmp/machine-id.XXXXXX)
        od -vAn -N16 -tx1 /dev/urandom | tr -d ' \\n' > "\$FAKE_MID"
        export APPTAINER_BINDPATH="\${FAKE_MID}:/etc/machine-id\${APPTAINER_BINDPATH:+,\$APPTAINER_BINDPATH}"
        export SINGULARITY_BINDPATH="\${FAKE_MID}:/etc/machine-id\${SINGULARITY_BINDPATH:+,\$SINGULARITY_BINDPATH}"
    fi

    schrodinger prime_mmgbsa \
        -job_type ${job_type} \
        -csv_output yes \
        -out_type PV \
        -jobname ${meta.id}_mmgbsa \
        -HOST localhost:${task.cpus} \
        -NJOBS ${task.cpus} \
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
  fi    fi    """
    stub:
    """
    touch ${meta.id}_mmgbsa-out.csv
    touch ${meta.id}_mmgbsa-out.maegz
    """
}