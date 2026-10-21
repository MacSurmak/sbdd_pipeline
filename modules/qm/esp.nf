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
    mult=\$(python3 -c "import json; print(json.load(open('${metadata_json}'))['spin_multiplicity'])")
    python3 ${projectDir}/bin/orca_esp.py \\
        --sdf ${refined_sdf} \\
        --charge ${meta.charge} \\
        --mult \${mult} \\
        --theory "${params.orca_esp_theory}" \\
        --basis "${params.orca_esp_basis}" \\
        --density "${params.orca_esp_density}" \\
        --scf_convergence "${params.orca_esp_scf_conv}" \\
        --scf_maxiter ${params.orca_esp_scf_maxiter} \\
        --nprocs ${task.cpus} \\
        --maxcore_mb ${params.orca_maxcore} \\
        --out_log "orca_esp.out"
    """
}
