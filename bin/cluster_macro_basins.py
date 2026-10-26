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
    parser = argparse.ArgumentParser(description="Partition trajectory phase space into macro-basins and identify medoids")
    parser.add_argument("--umap", required=True, help="Input umap_2d.npy array")
    parser.add_argument("--projections", required=True, help="Input projections.npy array")
    parser.add_argument("--coords", required=True, help="Input aligned_coords.npy array")
    parser.add_argument("--time_map", required=True, help="Input frame_time_map.json mapping")
    parser.add_argument(
        "--algorithm",
        default="hdbscan",
        choices=["hdbscan", "hierarchical", "auto_kmeans", "kmeans"],
        help="Clustering algorithm",
    )
    parser.add_argument("--n_clusters", type=int, default=3, help="Number of clusters (for kmeans)")
    parser.add_argument("--rmsd_cutoff", type=float, default=1.5, help="RMSD cutoff in Ångströms (for hierarchical)")
    parser.add_argument("--min_cluster_fraction", type=float, default=0.05, help="Minimum population fraction for a basin")
    parser.add_argument("--out_json", required=True, help="Output medoids_summary.json deliverable")
    parser.add_argument("--out_labels", required=True, help="Output cluster_labels.npy array")
    args = parser.parse_args()

    umap_p = Path(args.umap).resolve()
    proj_p = Path(args.projections).resolve()
    coords_p = Path(args.coords).resolve()
    time_p = Path(args.time_map).resolve()

    if not umap_p.exists() or not proj_p.exists() or not coords_p.exists() or not time_p.exists():
        raise FileNotFoundError("FATAL: Required input array(s) or time mapping missing for cluster_macro_basins!")

    umap_2d = np.load(umap_p)
    proj = np.load(proj_p)
    coords = np.load(coords_p)
    with open(time_p, "r", encoding="utf-8") as f:
        time_map = json.load(f)

    n_frames = umap_2d.shape[0]
    if n_frames != proj.shape[0] or n_frames != coords.shape[0] or n_frames != len(time_map):
        raise ValueError(
            f"FATAL: Inconsistent sample sizes! UMAP has {n_frames}, Projections has {proj.shape[0]}, "
            f"Coords has {coords.shape[0]}, Time map has {len(time_map)}"
        )

    reps_list = [item["replicate"] for item in time_map]
    time_list = [item["time_ns"] for item in time_map]
    reps_arr = np.array(reps_list)
    unique_reps = sorted(list(set(reps_list)))
    n_heavy_atoms = coords.shape[1] // 3

    algo = args.algorithm.lower()
    min_size = max(5, int(round(n_frames * float(args.min_cluster_fraction))))

    print(
        f"[CLUSTER_BASINS] Executing algorithm '{algo}' on {n_frames} frames "
        f"({n_heavy_atoms} heavy atoms, min_basin_size={min_size} frames)...",
        flush=True,
    )

    if algo == "hdbscan":
        try:
            from sklearn.cluster import HDBSCAN
        except ImportError:
            try:
                import hdbscan as hdb_mod
                HDBSCAN = hdb_mod.HDBSCAN
            except ImportError:
                raise ImportError(
                    "FATAL: Neither 'sklearn.cluster.HDBSCAN' nor 'hdbscan' module is available in environment!"
                )

        clusterer = HDBSCAN(min_cluster_size=min_size, min_samples=min_size, metric="euclidean")
        labels = clusterer.fit_predict(umap_2d)

        # Fallback if density was completely uniform and only noise was returned
        valid_cids = sorted([c for c in set(labels) if c >= 0])
        if not valid_cids:
            print("[CLUSTER_BASINS] Warning: HDBSCAN found 0 distinct clusters. Falling back to auto_kmeans...", flush=True)
            algo = "auto_kmeans"

    if algo == "hierarchical":
        from sklearn.cluster import AgglomerativeClustering

        # Exact pairwise Cartesian RMSD in Ångströms
        dists = cdist(coords, coords, metric="euclidean") / np.sqrt(float(n_heavy_atoms))
        agg = AgglomerativeClustering(
            metric="precomputed",
            linkage="average",
            distance_threshold=float(args.rmsd_cutoff),
            n_clusters=None,
        )
        labels = agg.fit_predict(dists)

    if algo == "auto_kmeans":
        from sklearn.metrics import silhouette_score

        max_k = min(8, max(2, n_frames // min_size))
        best_k = 2
        best_score = -1.0

        for k_cand in range(2, max_k + 1):
            km = KMeans(n_clusters=k_cand, n_init=10, random_state=42)
            c_labels = km.fit_predict(umap_2d)
            score = silhouette_score(umap_2d, c_labels)
            if score > best_score:
                best_score = score
                best_k = k_cand

        print(f"[CLUSTER_BASINS] Auto-K evaluated [2..{max_k}]: selected optimal K={best_k} (Silhouette={best_score:.3f})", flush=True)
        kmeans = KMeans(n_clusters=best_k, n_init=30, random_state=42)
        labels = kmeans.fit_predict(umap_2d)

    elif algo == "kmeans":
        k = max(1, min(args.n_clusters, n_frames))
        kmeans = KMeans(n_clusters=k, n_init=30, random_state=42)
        labels = kmeans.fit_predict(umap_2d)

    reps_list = [item["replicate"] for item in time_map]
    time_list = [item["time_ns"] for item in time_map]
    reps_arr = np.array(reps_list)
    unique_reps = sorted(list(set(reps_list)))

    raw_unique_cids = sorted([c for c in set(labels) if c >= 0])
    if not raw_unique_cids:
        raw_unique_cids = [0]
        labels = np.zeros(n_frames, dtype=np.int32)

    # First pass: identify initial medoids for genuine basins
    initial_stats = []
    for cid in raw_unique_cids:
        mask = (labels == cid)
        pop_count = int(np.sum(mask))
        if pop_count == 0:
            continue
        orig_indices = np.where(mask)[0]
        coords_sub = coords[mask]

        dists_rmsd = cdist(coords_sub, coords_sub, metric="euclidean") / np.sqrt(float(n_heavy_atoms))
        best_local_idx = np.argmin(np.sum(dists_rmsd, axis=1))
        medoid_frame = int(orig_indices[best_local_idx])
        spread = float(np.mean(dists_rmsd[best_local_idx]))
        initial_stats.append({
            "cid": cid,
            "medoid_frame": medoid_frame,
            "spread": spread,
        })

    # Zero Discard: Reassign noise frames (-1) to the nearest medoid by Cartesian heavy-atom RMSD
    noise_mask = (labels == -1)
    if np.any(noise_mask) and len(initial_stats) > 0:
        noise_indices = np.where(noise_mask)[0]
        medoid_frames = [s["medoid_frame"] for s in initial_stats]
        medoid_coords = coords[medoid_frames]
        for n_idx in noise_indices:
            pt = coords[n_idx : n_idx + 1]
            dists_to_medoids = cdist(pt, medoid_coords, metric="euclidean")[0] / np.sqrt(float(n_heavy_atoms))
            best_basin_idx = int(np.argmin(dists_to_medoids))
            labels[n_idx] = initial_stats[best_basin_idx]["cid"]

    # Final pass: exact 100% population accounting across all non-noise basins
    cluster_stats = []
    final_cids = sorted(list(set(labels)))
    for cid in final_cids:
        mask = (labels == cid)
        pop_count = int(np.sum(mask))
        if pop_count == 0:
            continue
        pop_pct = round((pop_count / n_frames) * 100.0, 2)
        orig_indices = np.where(mask)[0]

        rep_breakdown = {}
        for r in unique_reps:
            rep_breakdown[r] = int(np.sum((reps_arr == r) & mask))

        coords_sub = coords[mask]
        dists_rmsd = cdist(coords_sub, coords_sub, metric="euclidean") / np.sqrt(float(n_heavy_atoms))
        best_local_idx = np.argmin(np.sum(dists_rmsd, axis=1))
        medoid_frame = int(orig_indices[best_local_idx])
        spread = float(np.mean(dists_rmsd[best_local_idx]))

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
