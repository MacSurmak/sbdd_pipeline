/*
 * modules/schrodinger/split_prepared.nf
 * Splits fully prepared/protonated complex into separate receptor and ligand files.
 */

process SPLIT_PREPARED {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), path(complex_prepared_mae), path(poses_spec_json)

    output:
    tuple val(meta), path("${meta.name}_receptor.mae"), emit: prepared_mae
    tuple val(meta), path("${meta.name}_receptor.pdb"), emit: prepared_pdb
    tuple val(meta), path("*_pose.mae"), emit: pose_maes, optional: true
    tuple val(meta), path("*_pose.sdf"), emit: pose_sdfs, optional: true
    tuple val(meta), path("extracted_poses.json"), emit: poses_manifest

    script:
    """
    schrodinger run ${projectDir}/bin/split_prepared_complex.py \\
        --complex_mae ${complex_prepared_mae} \\
        --poses_spec_json ${poses_spec_json} \\
        --out_receptor_mae "${meta.name}_receptor.mae" \\
        --out_receptor_pdb "${meta.name}_receptor.pdb" \\
        --out_manifest extracted_poses.json
    """
}
