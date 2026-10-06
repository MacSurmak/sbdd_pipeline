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
    tuple val(meta), path("receptor.mae"), emit: prepared_mae
    tuple val(meta), path("receptor.pdb"), emit: prepared_pdb

    script:
    def lig_arg = meta.native_ligand ? "--native_ligand ${meta.native_ligand}" : ""
    """
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py light_prep \\
        --input_pdb ${input_pdb} \\
        --out_mae receptor.mae \\
        --out_pdb receptor.pdb \\
        ${lig_arg}
    """

    stub:
    """
    touch ${meta.id}_prepared.mae
    touch ${meta.id}_prepared.pdb
    """}
