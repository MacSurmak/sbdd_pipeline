/*
 * modules/schrodinger/ligprep.nf
 * Ligand Preparation with Epik and ForceField
 */

process LIGPREP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(ligand_file)
    val(target_ph)
    val(ph_tolerance)
    val(ff_version)
    val(max_states)

    output:
    tuple val(meta), path("${meta.id}_prepared.maegz"), emit: prepared_maegz
    tuple val(meta), path("${meta.id}_prepared.sdf"), emit: prepared_sdf

    script:
    def in_flag = (ligand_file.name.endsWith('.mae') || ligand_file.name.endsWith('.maegz')) ? "-imae" : "-isd"
    """
    schrodinger ligprep \\
        ${in_flag} ${ligand_file} \\
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

