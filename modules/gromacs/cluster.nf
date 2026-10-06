/*
 * modules/gromacs/cluster.nf
 * GROMACS trajectory clustering and automatic PDB splitting
 */

process GMX_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(clean_xtcs, stageAs: 'trajs_?/*'), path(ref_pdb), val(cutoff), val(method), val(dt), val(cluster_group)

    output:
    tuple val(meta), path("${meta.id}_centroid_*.pdb"), path("${meta.id}_centroid_*.mae"), emit: cluster_frames
    tuple val(meta), path("${meta.id}_pose_*.sdf"), path("${meta.id}_pose_*.mae"), emit: cluster_poses, optional: true
    path("${meta.id}_cluster_receipt.tsv"), emit: cluster_manifest

    script:
    def max_c = meta.max_clusters ?: 0
    def dt_arg = (dt && dt.toString() != '0') ? "--dt ${dt}" : ""
    """
    python3 ${projectDir}/bin/wrap_gromacs.py cluster \\
        --xtcs trajs_*/*.xtc \\
        --ref_pdb ${ref_pdb} \\
        --cutoff ${cutoff} \\
        --method ${method} \\
        --group ${cluster_group} \\
        --ligand "${meta.ligand ?: ''}" \\
        --max_clusters ${max_c} \\
        --out_id "${meta.id}" \\
        ${dt_arg}
    """

    stub:
    """
    touch ${meta.id}_centroid_1.pdb
    touch ${meta.id}_centroid_1.mae
    touch ${meta.id}_pose_1.sdf
    touch ${meta.id}_pose_1.mae
    printf "%s\\t%s/%s\\t%s/%s\\t%s/%s\\n" "${meta.id}_centroid_1" "\$PWD" "${meta.id}_centroid_1.pdb" "\$PWD" "${meta.id}_centroid_1.mae" "\$PWD" "${meta.id}_pose_1.sdf" > cluster_receipt_chunk.tsv
    """
}
