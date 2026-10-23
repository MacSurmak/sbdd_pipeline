/*
 * modules/analysis/cluster_basins.nf
 * Stage 2: K-Means macro-basin partitioning on 2D UMAP manifold
 * and medoid discovery in top-k projection space.
 */

process CLUSTER_BASINS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(umap_2d), path(projections), path(aligned_coords), path(frame_time_map)
    val(n_clusters)

    output:
    tuple val(meta), path("medoids_summary.json"), path("cluster_labels.npy"), emit: basins

    script:
    """
    python3 ${projectDir}/bin/cluster_macro_basins.py \\
        --umap ${umap_2d} \\
        --projections ${projections} \\
        --coords ${aligned_coords} \\
        --time_map ${frame_time_map} \\
        --n_clusters ${n_clusters} \\
        --out_json medoids_summary.json \\
        --out_labels cluster_labels.npy
    """
}
