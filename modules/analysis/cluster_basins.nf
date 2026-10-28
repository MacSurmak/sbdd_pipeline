/*
 * modules/analysis/cluster_basins.nf
 * Stage 2: K-Means macro-basin partitioning on 2D UMAP manifold
 * and medoid discovery in top-k projection space.
 */

process CLUSTER_BASINS {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(umap_2d), path(projections), path(aligned_coords), path(frame_time_map)
    val(algorithm)
    val(n_clusters)
    val(rmsd_cutoff)
    val(min_pop)

    output:
    tuple val(meta), path("medoids_summary.json"), path("cluster_labels.npy"), emit: basins

    script:
    """
    python3 ${projectDir}/bin/cluster_macro_basins.py \\
        --umap ${umap_2d} \\
        --projections ${projections} \\
        --coords ${aligned_coords} \\
        --time_map ${frame_time_map} \\
        --algorithm "${algorithm}" \\
        --n_clusters ${n_clusters} \\
        --rmsd_cutoff ${rmsd_cutoff} \\
        --min_cluster_pop_pct ${min_pop} \\
        --out_json medoids_summary.json \\
        --out_labels cluster_labels.npy
    """
}
