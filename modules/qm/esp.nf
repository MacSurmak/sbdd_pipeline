/*
 * modules/qm/esp.nf
 * Stage 2: Gas-phase canonical Merz-Kollman ESP single-point evaluation (HF/6-31G* TightSCF)
 */

process ORCA_ESP_CALC {
    tag "${meta.name}"

    input:
    tuple val(meta), path(refined_sdf), path(metadata_json)

    output:
    tuple val(meta), path(refined_sdf), path("orca_esp.out"), emit: esp_out

    script:
    """
    python3 ${projectDir}/bin/orca_esp.py \\
        --template "${projectDir}/assets/orca/esp.inp" \\
        --sdf ${refined_sdf} \\
        --charge ${meta.charge} \\
        --metadata_json ${metadata_json} \\
        --nprocs ${task.cpus} \\
        --maxcore_mb ${params.orca_maxcore} \\
        --out_log "orca_esp.out"
    """
}
