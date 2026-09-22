/*
 * modules/schrodinger/ligprep.nf
 * Ligand Preparation with Epik and ForceField
 */

process LIGPREP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(chunk_sdf)
    val(target_ph)
    val(ph_tolerance)
    val(ff_version)
    val(max_states)

    output:
    tuple val(meta), path("${meta.id}_prepared.maegz"), emit: prepared_maegz
    tuple val(meta), path("${meta.id}_prepared.sdf"), emit: prepared_sdf

    script:
    """
    schrodinger ligprep \\
        -isd ${chunk_sdf} \\
        -omae ${meta.id}_prepared.maegz \\
        -epik \\
        -ph ${target_ph} \\
        -pht ${ph_tolerance} \\
        -s ${max_states} \\
        -bff ${ff_version} \\
        -g \\
        -WAIT

    schrodinger structconvert ${meta.id}_prepared.maegz ${meta.id}_prepared.sdf
    """

    stub:
    """
    touch ${meta.id}_prepared.maegz
    touch ${meta.id}_prepared.sdf
    """
}

