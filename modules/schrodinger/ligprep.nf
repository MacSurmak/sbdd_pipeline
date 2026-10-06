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
    tuple val(meta), path("ligand.maegz"), emit: prepared_maegz
    tuple val(meta), path("ligand.sdf"), emit: prepared_sdf

    script:
    """
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py ligprep \\
        --ligand_file ${ligand_file} \\
        --target_ph ${target_ph} \\
        --ph_tolerance ${ph_tolerance} \\
        --ff_version ${ff_version} \\
        --max_states ${max_states} \\
        --out_maegz ligand.maegz \\
        --out_sdf ligand.sdf
    """

    stub:
    """
    touch ${meta.id}_prepared.maegz
    touch ${meta.id}_prepared.sdf
    """
}

