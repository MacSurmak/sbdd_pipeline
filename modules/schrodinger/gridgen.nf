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

    output:
    tuple val(meta), path("${meta.id}_grid.zip"), emit: grid_zip

    script:
    """
    cat << EOF > gridgen.in
GRIDFILE ${meta.id}_grid.zip
RECEP_FILE ${receptor_mae}
REF_LIGAND_FILE ${ref_lig_mae}
INNERBOX ${inner_box}
OUTERBOX ${outer_box}
FORCEFIELD ${forcefield}
EOF

    schrodinger glide gridgen.in -OVERWRITE -WAIT
    """

    stub:
    """
    touch ${meta.id}_grid.zip
    """
}

