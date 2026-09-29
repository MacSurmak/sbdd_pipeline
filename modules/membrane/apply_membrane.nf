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
    tuple val(meta), path("${meta.id}_prepared_mem.mae"), emit: prepared_mae

    script:
    """
    schrodinger run ${projectDir}/bin/apply_membrane.py \\
        --input ${prepared_mae} \\
        --info ${membrane_info} \\
        --output tmp_${meta.id}_prepared_mem.mae

    mv tmp_${meta.id}_prepared_mem.mae ${meta.id}_prepared_mem.mae
    """

    stub:
    """
    touch ${meta.id}_prepared_mem.mae
    """
}
