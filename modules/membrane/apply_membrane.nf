/*
 * modules/membrane/apply_membrane.nf
 * Implicit Membrane Property Injection for Prime
 */

process APPLY_MEMBRANE {
    tag "${meta.id}"
    label "process_single"
    label "process_master"

    input:
    tuple val(meta), path(prepared_mae), path(membrane_info), path(apply_script)

    output:
    tuple val(meta), path("${meta.id}_membrane.mae"), emit: membrane_mae

    script:
    """
    schrodinger run python3 ${apply_script} \\
        --input ${prepared_mae} \\
        --info ${membrane_info} \\
        --output ${meta.id}_membrane.mae
    """

    stub:
    """
    touch ${meta.id}_membrane.mae
    """
}
