/*
 * modules/schrodinger/split_prepared.nf
 * Splits fully prepared/protonated complex into separate receptor and ligand files.
 */

process SPLIT_PREPARED {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), path(complex_prepared_mae)

    output:
    tuple val(meta), path("${meta.id}_prepared.mae"), emit: prepared_mae
    tuple val(meta), path("${meta.id}_prepared.pdb"), emit: prepared_pdb
    tuple val(meta), path("${meta.id}_ref_native_prepared.mae"), emit: ref_native_mae
    tuple val(meta), path("${meta.id}_ref_native_prepared.sdf"), emit: ref_native_sdf

    script:
    def lig_arg = meta.native_ligand ? "--ligand ${meta.native_ligand}" : ""
    """
    schrodinger run ${projectDir}/bin/split_prepared_complex.py \\
        --complex_mae ${complex_prepared_mae} \\
        ${lig_arg} \\
        --out_receptor_mae ${meta.id}_prepared.mae \\
        --out_receptor_pdb ${meta.id}_prepared.pdb \\
        --out_ligand_mae ${meta.id}_ref_native_prepared.mae \\
        --out_ligand_sdf ${meta.id}_ref_native_prepared.sdf
    """

    stub:
    """
    touch ${meta.id}_prepared.mae
    touch ${meta.id}_prepared.pdb
    touch ${meta.id}_ref_native_prepared.mae
    touch ${meta.id}_ref_native_prepared.sdf
    """
}
