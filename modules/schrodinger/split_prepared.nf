/*
 * modules/schrodinger/split_prepared.nf
 * Splits fully prepared/protonated complex into separate receptor and ligand files.
 */

process SPLIT_PREPARED {
    tag "${meta.id}"
    label "process_single"
    publishDir path: { "${params.outdir}/prepared/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(complex_prepared_mae)

    output:
    tuple val(meta), path("receptor.mae"), emit: prepared_mae
    tuple val(meta), path("receptor.pdb"), emit: prepared_pdb
    tuple val(meta), path("native.mae"), emit: ref_native_mae, optional: true
    tuple val(meta), path("native.sdf"), emit: ref_native_sdf, optional: true

    output:
    tuple val(meta), path("receptor.mae"), emit: prepared_mae
    tuple val(meta), path("receptor.pdb"), emit: prepared_pdb
    tuple val(meta), path("native.mae"), emit: ref_native_mae, optional: true
    tuple val(meta), path("native.sdf"), emit: ref_native_sdf, optional: true
    tuple val(meta), path("pose_*.sdf"), emit: cofactor_sdfs, optional: true
    tuple val(meta), path("pose_*.mae"), emit: cofactor_maes, optional: true
    tuple val(meta), path("extracted_poses.json"), emit: poses_manifest

    script:
    def lig_arg = meta.native_ligand ? "--ligand ${meta.native_ligand}" : ""
    def cof_list = meta.cofactors ? (meta.cofactors instanceof List ? meta.cofactors.join(' ') : meta.cofactors) : ""
    def cof_arg = cof_list ? "--cofactors ${cof_list}" : ""
    """
    schrodinger run ${projectDir}/bin/split_prepared_complex.py \\
        --complex_mae ${complex_prepared_mae} \\
        ${lig_arg} \\
        ${cof_arg} \\
        --out_receptor_mae receptor.mae \\
        --out_receptor_pdb receptor.pdb \\
        --out_ligand_mae native.mae \\
        --out_ligand_sdf native.sdf \\
        --out_manifest extracted_poses.json
    """

    stub:
    """
    touch receptor.mae
    touch receptor.pdb
    touch native.mae
    touch native.sdf
    echo "[]" > extracted_poses.json
    """
}
