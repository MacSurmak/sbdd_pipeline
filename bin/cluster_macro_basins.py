#!/usr/bin/env python3
"""
bin/cluster_macro_basins.py
Atomic Stage 2: K-Means macro-basin partitioning in 2D UMAP manifold space
and physical medoid discovery in top-k PCA/TICA subspace (0% artificial noise discard).
"""

import argparse
import json
import sys
from pathlib import Path
import numpy as np
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist


def main():
    parser = argparse.ArgumentParser(description="Partition 2D manifold into macro-basins and identify medoids")
    parser.add_argument("--umap", required=True, help="Input umap_2d.npy array")
    parser.add_argument("--projections", required=True, help="Input projections.npy array")
    parser.add_argument("--time_map", required=True, help="Input frame_time_map.json mapping")
    parser.add_argument("--n_clusters", type=int, required=True, help="Number of macro-basins (K)")
    parser.add_argument("--out_json", required=True, help="Output medoids_summary.json deliverable")
    parser.add_argument("--out_labels", required=True, help="Output cluster_labels.npy array")
    args = parser.parse_args()

    umap_p = Path(args.umap).resolve()
    proj_p = Path(args.projections).resolve()
    time_p = Path(args.time_map).resolve()

    if not umap_p.exists() or not proj_p.exists() or not time_p.exists():
        raise FileNotFoundError("FATAL: Required input array(s) or time mapping missing for cluster_macro_basins!")

    umap_2d = np.load(umap_p)
    proj = np.load(proj_p)
    with open(time_p, "r", encoding="utf-8") as f:
        time_map = json.load(f)

    n_frames = umap_2d.shape[0]
    if n_frames != proj.shape[0] or n_frames != len(time_map):
        raise ValueError(
            f"FATAL: Inconsistent sample sizes! UMAP has {n_frames}, Projections has {proj.shape[0]}, "
            f"Time map has {len(time_map)}"
        )

    k = max(1, min(args.n_clusters, n_frames))
    print(f"[CLUSTER_BASINS] Partitioning {n_frames} frames into K={k} physical macro-basins...", flush=True)

    kmeans = KMeans(n_clusters=k, n_init=30, random_state=42)
    labels = kmeans.fit_predict(umap_2d)

    reps_list = [item["replicate"] for item in time_map]
    time_list = [item["time_ns"] for item in time_map]
    reps_arr = np.array(reps_list)
    unique_reps = sorted(list(set(reps_list)))

    cluster_stats = []
    for cid in range(k):
        mask = labels == cid
        pop_count = int(np.sum(mask))
        if pop_count == 0:
            continue
        pop_pct = round((pop_count / n_frames) * 100.0, 2)
        orig_indices = np.where(mask)[0]

        rep_breakdown = {}
        for r in unique_reps:
            rep_breakdown[r] = int(np.sum((reps_arr == r) & mask))

        # True medoid: frame minimizing sum of Euclidean distances to all other basin members in PCA space
        proj_sub = proj[mask]
        dists = cdist(proj_sub, proj_sub, metric="euclidean")
        best_local_idx = np.argmin(np.sum(dists, axis=1))
        medoid_frame = int(orig_indices[best_local_idx])
        spread = float(np.mean(dists[best_local_idx]))

        cluster_stats.append({
            "raw_cluster_id": int(cid),
            "population_count": pop_count,
            "population_pct": pop_pct,
            "medoid_frame": medoid_frame,
            "medoid_time_ns": time_list[medoid_frame],
            "medoid_replicate": reps_list[medoid_frame],
            "compactness_spread": round(spread, 3),
            "replicate_breakdown": rep_breakdown,
            "umap_coords": [round(float(umap_2d[medoid_frame, 0]), 2), round(float(umap_2d[medoid_frame, 1]), 2)],
            "pca_coords": [round(float(proj[medoid_frame, 0]), 2), round(float(proj[medoid_frame, 1]), 2)],
        })

    cluster_stats.sort(key=lambda x: x["population_pct"], reverse=True)
    for rank, stat in enumerate(cluster_stats, start=1):
        stat["rank"] = rank

    np.save(Path(args.out_labels).resolve(), labels)
    with open(Path(args.out_json).resolve(), "w", encoding="utf-8") as f:
        json.dump(cluster_stats, f, indent=2)

    print(f"✓ Macro-basin clustering completed: {len(cluster_stats)} basins identified (0% noise).", flush=True)


if __name__ == "__main__":
    main()
