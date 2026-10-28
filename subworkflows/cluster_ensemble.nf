/*
 * subworkflows/cluster_ensemble.nf
 * Pure analytical phase-space clustering (PCA/TICA + UMAP -> Basins -> Reporting)
 * Decoupled from physical structure extraction: 100% focused on coordinate partitioning.
 */

include { TRAJ_DIM_REDUCTION } from '../modules/analysis/dim_reduction'
include { CLUSTER_BASINS     } from '../modules/analysis/cluster_basins'
include { RENDER_REPORTS     } from '../modules/analysis/render_reports'

workflow CLUSTER_ENSEMBLE {
    take:
    ch_cluster_input // tuple val(meta), path(ref_pdb), path(xtc_files), path(ndx_file)
    align_group
    cluster_group
    dim_reduction
    algorithm
    start_ns
    sample_interval_ps
    tica_lag_ns
    variance_cutoff
    pocket_radius
    n_clusters
    rmsd_cutoff
    min_pop
    temp_k

    main:
    // 1. Unified Slicing Stage: Stream and downsample trajectories synchronized with MM-PBSA
    TRAJ_DIM_REDUCTION(
        ch_cluster_input,
        align_group,
        cluster_group,
        dim_reduction,
        start_ns,
        sample_interval_ps,
        tica_lag_ns,
        variance_cutoff,
        pocket_radius
    )

    // 2. Fast Stage: Partition into macro-basins without guessing K
    ch_basins_in = TRAJ_DIM_REDUCTION.out.reduced_data.map { meta, proj, umap, eig, coords, time_map, feat_meta ->
        [ meta, umap, proj, coords, time_map ]
    }
    CLUSTER_BASINS(ch_basins_in, algorithm, n_clusters, rmsd_cutoff, min_pop)

    // 3. Reporting Stage: FES computation, Plotly HTML, publication PNG, PyMOL session
    ch_render_in = ch_cluster_input.map { meta, ref_pdb, xtcs, ndx -> [ meta, ref_pdb ] }
        .join(TRAJ_DIM_REDUCTION.out.reduced_data.map { meta, proj, umap, eig, coords, time_map, feat_meta -> [ meta, umap, proj, time_map, feat_meta ] }, by: 0)
        .join(CLUSTER_BASINS.out.basins, by: 0)
        .map { meta, ref_pdb, umap, proj, time_map, feat_meta, medoids_json, labels ->
            [ meta, ref_pdb, umap, proj, labels, medoids_json, time_map, feat_meta ]
        }
    RENDER_REPORTS(ch_render_in, align_group, cluster_group, dim_reduction, algorithm, temp_k)

    emit:
    manifest = RENDER_REPORTS.out.reports.map { meta, html, png, pml, ref, manifest_json ->
        [ meta, manifest_json ]
    }
}
