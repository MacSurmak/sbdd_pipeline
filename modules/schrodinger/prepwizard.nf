/*
 * modules/schrodinger/prepwizard.nf
 * Protein Preparation Wizard
 */

process PREPWIZARD {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(complex_mae)
    val(target_ph)
    val(ph_tolerance)
    val(min_rmsd)
    val(cap_termini)

    output:
    tuple val(meta), path("receptor.mae"), emit: prepared_mae

    script:
    def cap_flag = cap_termini ? "-captermini" : ""
    """
    schrodinger prepwizard \\
        -fix -fillsidechains \\
        -propka_pH ${target_ph} \\
        -epik_pH ${target_ph} \\
        -epik_pHt ${ph_tolerance} \\
        -rmsd ${min_rmsd} \\
        -preserve_st_titles \\
        ${cap_flag} \\
        -WAIT \\
        ${complex_mae} \\
        receptor.mae
    """

    stub:
    """
    touch ${meta.id}_prepared.mae
    """
}

