/*
 * modules/analysis/dim_reduction.nf
 * Stage 1: In-memory coordinate streaming, on-the-fly Kabsch superimposition,
 * and dimensionality reduction (PCA/TICA + 2D UMAP).
 */

process TRAJ_DIM_REDUCTION {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(ref_pdb), path(xtc_files, stageAs: 'trajs/*'), path(ndx_file)
    val(align_group)
    val(cluster_group)
    val(dim_reduction)
    val(tica_lag_ns)
    val(variance_cutoff)
    val(pocket_radius)

    output:
    tuple val(meta), path("projections.npy"), path("umap_2d.npy"), path("eigenvalues.npy"), path("frame_time_map.json"), path("feature_meta.json"), emit: reduced_data

    script:
    def lig_arg = meta.ligand ? "--ligand_resname ${meta.ligand}" : ""
    def ndx_arg = (ndx_file && ndx_file.name != 'NO_FILE') ? "--ndx ${ndx_file}" : ""
    """
    python3 ${projectDir}/bin/traj_dim_reduction.py \\
        --struct ${ref_pdb} \\
        --xtc trajs/*.xtc \\
        ${ndx_arg} \\
        --align_group "${align_group}" \\
        --cluster_group "${cluster_group}" \\
        ${lig_arg} \\
        --pocket_radius ${pocket_radius} \\
        --dim_reduction ${dim_reduction} \\
        --tica_lag_ns ${tica_lag_ns} \\
        --variance_cutoff ${variance_cutoff} \\
        --outdir .
    """
}
