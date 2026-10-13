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
    def cap_flag = cap_termini ? "--cap_termini" : ""
    """
    python3 ${projectDir}/bin/wrap_schrodinger.py prepwizard \\
        --complex_mae ${complex_mae} \\
        --out_mae receptor.mae \\
        --target_ph ${target_ph} \\
        --ph_tolerance ${ph_tolerance} \\
        --min_rmsd ${min_rmsd} \\
        ${cap_flag}
    """

    stub:
    """
    touch ${meta.id}_prepared.mae
    """
}

