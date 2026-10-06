/*
 * modules/schrodinger/gridgen.nf
 * Receptor Grid Generation Module for Glide Docking
 */

process GLIDE_GRIDGEN {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(receptor_mae), path(ref_lig_mae)
    val(inner_box)
    val(outer_box)
    val(forcefield)
    val(recep_vscale)

    output:
    tuple val(meta), path("${meta.id}_grid.zip"), emit: grid_zip

    script:
    """
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py gridgen \\
        --receptor_mae ${receptor_mae} \\
        --ref_lig_mae ${ref_lig_mae} \\
        --inner_box "${inner_box}" \\
        --outer_box "${outer_box}" \\
        --forcefield ${forcefield} \\
        --recep_vscale ${recep_vscale} \\
        --out_zip ${meta.id}_grid.zip
    """

    stub:
    """
    touch ${meta.id}_grid.zip
    """
}
