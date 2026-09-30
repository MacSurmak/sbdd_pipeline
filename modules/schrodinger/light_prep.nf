/*
 * modules/schrodinger/light_prep.nf
 * Light preparation for MD frames/clusters.
 * Reassigns OPLS4 charges without altering heavy atom coordinates or protonation states.
 */

process LIGHT_PREPWIZARD {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(input_pdb)

    output:
    tuple val(meta), path("${meta.id}_lightprep.mae"), emit: prepared_mae

    script:
    """
    schrodinger prepwizard \\
        -reassign_charges -noimpref -nopka -delwater_nocons \\
        -WAIT ${input_pdb} ${meta.id}_lightprep.mae
    """

    stub:
    """
    touch ${meta.id}_lightprep.mae
    """
}
