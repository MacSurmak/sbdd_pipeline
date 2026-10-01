/*
 * modules/schrodinger/gridgen.nf
 * Receptor Grid Generation Module for Glide Docking
 */

process GLIDE_GRIDGEN {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(receptor_mae, stageAs: 'rec_dir/*'), path(ref_lig_mae, stageAs: 'ref_dir/*')
    val(inner_box)
    val(outer_box)
    val(forcefield)

    output:
    tuple val(meta), path("${meta.id}_grid.zip"), emit: grid_zip

    script:
    """
    REF_BASE=\$(basename "${ref_lig_mae}")
    REF_MAE="ref_ligand.mae"
    if [[ "\$REF_BASE" == *.mae || "\$REF_BASE" == *.maegz ]]; then
        REF_MAE="${ref_lig_mae}"
    elif [[ "\$REF_BASE" != "NO_FILE" && -f "${ref_lig_mae}" ]]; then
        schrodinger structconvert "${ref_lig_mae}" ref_ligand.mae
        REF_MAE="ref_ligand.mae"
    else
        REF_MAE="${ref_lig_mae}"
    fi

    cat << EOF > gridgen.in
GRIDFILE ${meta.id}_grid.zip
RECEP_FILE ${receptor_mae}
REF_LIGAND_FILE \${REF_MAE}
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

