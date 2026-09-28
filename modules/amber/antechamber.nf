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
    def name = meta.mol_id
    """
    charge=\$(cat ${charge_file} | tr -d '[:space:]')
    
    echo "Running Antechamber for ${name} with net charge \$charge"
    antechamber \
        -i ${mol2} \
        -fi mol2 \
        -o ${name}_gaff2.mol2 \
        -fo mol2 \
        -at gaff2 \
        -c bcc \
        -nc "\$charge" \
        -rn ${name}

    antechamber -i ${name}_gaff2.mol2 -fi mol2 -o ${name}_gaff2.pdb -fo pdb

    parmchk2 \
        -i ${name}_gaff2.mol2 \
        -f mol2 \
        -o ${name}.frcmod \
        -s 2
    """

    stub:
    """
    touch UNK_gaff2.mol2
    touch UNK.frcmod
    """
}

