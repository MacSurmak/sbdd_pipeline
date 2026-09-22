/*
 * modules/membrane/apply_membrane.nf
 * Implicit Membrane Property Injection for Prime
 * Self-contained worker process: zero external script file staging.
 */

process APPLY_MEMBRANE {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), path(prepared_mae), path(membrane_info)

    output:
    tuple val(meta), path("${meta.id}_membrane.mae"), emit: membrane_mae

    script:
    """
    schrodinger run python3 ${projectDir}/bin/apply_membrane.py \\
        --input ${prepared_mae} \\
        --info ${membrane_info} \\
        --output ${meta.id}_membrane.mae
    """

    stub:
    """
    touch ${meta.id}_membrane.mae
    """
}
