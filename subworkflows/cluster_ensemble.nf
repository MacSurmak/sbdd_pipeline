/*
 * subworkflows/cluster_ensemble.nf
 * Orchestrates 4 atomic Slurm tasks for conformational landscape clustering
 */

include { TRAJ_DIM_REDUCTION } from '../modules/analysis/dim_reduction'
include { CLUSTER_BASINS     } from '../modules/analysis/cluster_basins'
include { EXTRACT_CENTROIDS  } from '../modules/analysis/extract_centroids'
include { RENDER_REPORTS     } from '../modules/analysis/render_reports'

workflow CLUSTER_ENSEMBLE {
    take:
    ch_cluster_input // tuple val(meta), path(ref_pdb), path(xtc_files), path(ndx_file), path(template_sdf)
    align_group
    cluster_group
    dim_reduction
    tica_lag_ns
    variance_cutoff
    pocket_radius
    n_clusters
    temp_k

    main:
    // 1. Heavy Stage: Stream trajectories and project onto PCA/TICA + UMAP (Cached by Nextflow)
    ch_dim_in = ch_cluster_input.map { meta, ref_pdb, xtcs, ndx, tmpl_sdf ->
        [ meta, ref_pdb, xtcs, ndx ]
    }
    TRAJ_DIM_REDUCTION(
        ch_dim_in,
        ch_cluster_input.map { meta, ref_pdb, xtcs, ndx, tmpl_sdf -> meta.replicate_names },
        align_group,
        cluster_group,
        dim_reduction,
        tica_lag_ns,
        variance_cutoff,
        pocket_radius
    )

    // 2. Fast Stage: Partition into K macro-basins in UMAP manifold
    ch_basins_in = TRAJ_DIM_REDUCTION.out.reduced_data.map { meta, proj, umap, eig, time_map, feat_meta ->
        [ meta, umap, proj, time_map ]
    }
    CLUSTER_BASINS(ch_basins_in, n_clusters)

    // 3. Coordinate Stage: Direct seek frame extraction of physical medoids
    ch_extract_in = ch_cluster_input
        .join(CLUSTER_BASINS.out.basins, by: 0)
        .map { meta, ref_pdb, xtcs, ndx, tmpl_sdf, medoids_json, labels ->
            [ meta, ref_pdb, xtcs, medoids_json, tmpl_sdf ]
        }
    EXTRACT_CENTROIDS(ch_extract_in)

    // 4. Reporting Stage: FES computation, Plotly HTML, PNG, and PyMOL inspection script
    ch_render_in = ch_cluster_input.map { meta, ref_pdb, xtcs, ndx, tmpl_sdf -> [ meta, ref_pdb ] }
        .join(TRAJ_DIM_REDUCTION.out.reduced_data.map { meta, proj, umap, eig, time_map, feat_meta -> [ meta, umap, proj, time_map ] }, by: 0)
        .join(CLUSTER_BASINS.out.basins, by: 0)
        .map { meta, ref_pdb, umap, proj, time_map, medoids_json, labels ->
            [ meta, ref_pdb, umap, proj, labels, medoids_json, time_map ]
        }
    RENDER_REPORTS(ch_render_in, temp_k)

    emit:
    manifest = RENDER_REPORTS.out.reports.map { meta, html, png, pml, ref, manifest_json ->
        [ meta, manifest_json ]
    }
}
