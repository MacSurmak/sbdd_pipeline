/*
 * modules/amber/antechamber.nf
 * Antechamber AM1-BCC Charge Calculation and GAFF2 Parametrization
 */

process ANTECHAMBER {
    tag "${meta.id}_${meta.mol_id}"
    label 'process_single'

    input:
    tuple val(meta), path(mol2), path(charge_file)

    output:
    tuple val(meta), path("${meta.mol_id}_gaff2.mol2"), path("${meta.mol_id}.frcmod"), path("${meta.mol_id}_gaff2.pdb"), emit: gaff_params

    script:
    """
    python3 ${projectDir}/bin/wrap_amber.py antechamber \\
        --mol2 ${mol2} \\
        --charge_file ${charge_file} \\
        --mol_id ${meta.mol_id}
    """

    stub:
    """
    touch UNK_gaff2.mol2
    touch UNK.frcmod
    """
}

