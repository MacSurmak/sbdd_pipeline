/*
 * modules/analysis/render_reports.nf
 * Stage 4: Boltzmann Free Energy Surface calculation, interactive Plotly dashboard,
 * publication PNG, and portable relative PyMOL session.
 */

process RENDER_REPORTS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(ref_pdb), path(umap_2d), path(projections), path(cluster_labels), path(medoids_summary), path(frame_time_map), path(feature_meta)
    val(align_group)
    val(cluster_group)
    val(dim_reduction)
    val(temp_k)

    output:
    tuple val(meta), path("macro_clustering_interactive.html"), path("macro_clustering_landscape.png"), path("inspect_clusters.pml"), path("reference.pdb"), path("clusters_manifest.json"), emit: reports

    script:
    """
    python3 ${projectDir}/bin/render_cluster_reports.py \\
        --struct ${ref_pdb} \\
        --umap ${umap_2d} \\
        --projections ${projections} \\
        --labels ${cluster_labels} \\
        --medoids_json ${medoids_summary} \\
        --time_map ${frame_time_map} \\
        --feature_meta ${feature_meta} \\
        --align_group "${align_group}" \\
        --cluster_group "${cluster_group}" \\
        --dim_reduction "${dim_reduction}" \\
        --temp_k ${temp_k} \\
        --outdir .
    """
}
