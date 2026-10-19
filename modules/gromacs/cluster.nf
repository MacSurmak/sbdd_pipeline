/*
 * modules/gromacs/cluster.nf
 * GROMACS trajectory clustering and automatic PDB splitting
 */

process GMX_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(clean_xtcs, stageAs: 'trajs_?/*'), path(ref_pdb), val(cutoff), val(method), val(dt), val(cluster_group)
    val(water_names)
    val(ion_names)
    val(lipid_names)

    output:
    tuple val(meta), path("*_receptor.pdb"), path("*_receptor.mae"), emit: cluster_frames
    tuple val(meta), path("*_view.pdb"), emit: view_pdbs, optional: true
    tuple val(meta), path("*_pose.sdf"), path("*_pose.mae"), emit: cluster_poses, optional: true
    path("${meta.id}_cluster_receipt.json"), emit: cluster_manifest

    script:
    def max_c = meta.max_clusters ?: 0
    def dt_arg = (dt && dt.toString() != '0') ? "--dt ${dt}" : ""
    def sys_identifier = meta.system ?: meta.id
    """
    python3 ${projectDir}/bin/wrap_gromacs.py cluster \\
        --xtcs trajs_*/*.xtc \\
        --ref_pdb ${ref_pdb} \\
        --cutoff ${cutoff} \\
        --method ${method} \\
        --group ${cluster_group} \\
        --ligand "${meta.ligand ?: ''}" \\
        --max_clusters ${max_c} \\
        --sys_name "${sys_identifier}" \\
        --water_names "${water_names}" \\
        --ion_names "${ion_names}" \\
        --lipid_names "${lipid_names}" \\
        --out_receipt "${meta.id}_cluster_receipt.json" \\
        ${dt_arg}
    """
}
