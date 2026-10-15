/*
 * modules/schrodinger/split_prepared.nf
 * Splits fully prepared/protonated complex into separate receptor and ligand files.
 */

process SPLIT_PREPARED {
    tag "${meta.id}"
    label "process_single"
    publishDir path: { "${params.outdir}/prepared_frames/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(complex_prepared_mae)

    output:
    tuple val(meta), path("receptor.mae"), emit: prepared_mae
    tuple val(meta), path("receptor.pdb"), emit: prepared_pdb
    tuple val(meta), path("native.mae"), emit: ref_native_mae, optional: true
    tuple val(meta), path("native.sdf"), emit: ref_native_sdf, optional: true

    script:
    def lig_arg = meta.native_ligand ? "--ligand ${meta.native_ligand}" : ""
    """
    schrodinger run ${projectDir}/bin/split_prepared_complex.py \\
        --complex_mae ${complex_prepared_mae} \\
        ${lig_arg} \\
        --out_receptor_mae receptor.mae \\
        --out_receptor_pdb receptor.pdb \\
        --out_ligand_mae native.mae \\
        --out_ligand_sdf native.sdf
    """

    stub:
    """
    touch ${meta.id}_prepared.mae
    touch ${meta.id}_prepared.pdb
    touch ${meta.id}_ref_native_prepared.mae
    touch ${meta.id}_ref_native_prepared.sdf
    """
}
