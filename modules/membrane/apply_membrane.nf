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
    tuple val(meta), path("${meta.id}_prepared.mae"), emit: prepared_mae
    tuple val(meta), path("${meta.id}_prepared.pdb"), emit: prepared_pdb

    script:
    """
    cp \$(which apply_membrane.py) .
    schrodinger run apply_membrane.py \\
        --input ${prepared_mae} \\
        --info ${membrane_info} \\
        --output tmp_${meta.id}_prepared.mae

    mv tmp_${meta.id}_prepared.mae ${meta.id}_prepared.mae
    schrodinger structconvert ${meta.id}_prepared.mae ${meta.id}_prepared.pdb
    """

    stub:
    """
    touch ${meta.id}_prepared.mae
    touch ${meta.id}_prepared.pdb
    """
}
