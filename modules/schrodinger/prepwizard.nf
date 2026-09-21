/*
 * modules/schrodinger/prepwizard.nf
 * Protein Preparation Wizard
 */

process PREPWIZARD {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_master'

    input:
    tuple val(meta), path(raw_rec)
    val(target_ph)
    val(ph_tolerance)
    val(min_rmsd)
    val(cap_termini)

    output:
    tuple val(meta), path("${meta.id}_prepared.mae"), emit: prepared_mae
    tuple val(meta), path("${meta.id}_prepared.pdb"), emit: prepared_pdb

    script:
    def cap_flag = cap_termini ? "-captermini" : ""
    """
    schrodinger prepwizard \
        -fix -fillsidechains \
        -propka_pH ${target_ph} \
        -epik_pH ${target_ph} \
        -epik_pHt ${ph_tolerance} \
        -rmsd ${min_rmsd} \
        -preserve_st_titles \
        ${cap_flag} \
        -WAIT \
        ${raw_rec} \
        ${meta.id}_prepared.mae

    schrodinger structconvert ${meta.id}_prepared.mae ${meta.id}_prepared.pdb
    """

    stub:
    """
    touch ${meta.id}_prepared.mae
    touch ${meta.id}_prepared.pdb
    """
}

