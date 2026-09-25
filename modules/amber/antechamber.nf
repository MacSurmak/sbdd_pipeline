/*
 * modules/amber/antechamber.nf
 * Antechamber AM1-BCC Charge Calculation and GAFF2 Parametrization
 */

process ANTECHAMBER {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(ligand_mol2), path(charge_file)

    output:
    tuple val(meta), path("UNK_gaff2.mol2"), path("UNK.frcmod"), emit: gaff_params

    script:
    """
    charge=\$(cat ${charge_file} | tr -d '[:space:]')
    [ -z "\$charge" ] && charge=0

    echo "Running Antechamber with net charge \$charge"
    antechamber \
        -i ${ligand_mol2} \
        -fi mol2 \
        -o UNK_gaff2.mol2 \
        -fo mol2 \
        -at gaff2 \
        -c bcc \
        -nc "\$charge" \
        -rn UNK

    parmchk2 \
        -i UNK_gaff2.mol2 \
        -f mol2 \
        -o UNK.frcmod \
        -s 2
    """

    stub:
    """
    touch UNK_gaff2.mol2
    touch UNK.frcmod
    """
}

