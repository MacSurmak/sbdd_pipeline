#!/usr/bin/env python3
"""
bin/stratify_mmpbsa.py
Conformation-conditioned MM-PBSA stratification engine.
Performs Point-to-State projection between independent thermodynamic evaluations
Delta_G(t_j) and conformational cluster states C(t_j) without synthetic interpolation.
Computes Boltzmann ensemble aggregate and renders multi-replica color-coded RMSD trajectories.
"""

import argparse
import json
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def parse_args():
    parser = argparse.ArgumentParser(
        description="Stratify MM-PBSA binding free energy across conformational macro-basins"
    )
    parser.add_argument("--system", required=True, help="Canonical system name")
    parser.add_argument("--labels", required=True, help="Path to cluster_labels.npy")
    parser.add_argument("--time_map", required=True, help="Path to frame_time_map.json")
    parser.add_argument(
        "--mmpbsa_dir",
        required=True,
        help="Directory containing replicate MM-PBSA outputs",
    )
    parser.add_argument(
        "--metrics_dir",
        default="",
        help="Optional base directory to search for rmsd.csv files",
    )
    parser.add_argument(
        "--temp_k", type=float, default=310.15, help="Simulation temperature in Kelvin"
    )
    parser.add_argument(
        "--out_json", default="stratification_summary.json", help="Output summary JSON"
    )
    parser.add_argument(
        "--out_plot",
        default="stratified_rmsd_plot.png",
        help="Output multi-replica plot PNG",
    )
    return parser.parse_args()


def load_replicate_mmpbsa(
    mmpbsa_dir: Path, system_name: str
) -> dict[str, pd.DataFrame]:
    """Loads and sanitizes per-replicate DELTA binding free energy time series."""
    rep_dfs = {}
    target_dirs = sorted(
        [d for d in mmpbsa_dir.glob("*") if system_name in d.name and d.is_dir()]
    )
    if not target_dirs:
        raise FileNotFoundError(
            f"FATAL: No replicate output folders found for system '{system_name}' in {mmpbsa_dir}!"
        )

    for r_dir in target_dirs:
        csv_f = r_dir / "thermodynamics_mmpbsa_pb.csv"
        if not csv_f.exists() or csv_f.stat().st_size == 0:
            continue

        df = pd.read_csv(csv_f)
        target_col = next((c for c in df.columns if "TOTAL" in c.upper()), None)
        if not target_col:
            raise KeyError(f"FATAL: No 'TOTAL' energy column found in {csv_f}!")

        # Clean frame index column
        frame_col = next((c for c in df.columns if "FRAME" in c.upper()), df.columns[0])

        df_clean = pd.DataFrame(
            {
                "sample_idx": np.arange(len(df)),
                "energy_frame": pd.to_numeric(df[frame_col], errors="coerce")
                .fillna(0)
                .astype(int),
                "delta_g": pd.to_numeric(df[target_col], errors="coerce"),
            }
        ).dropna()

        # Deduce replicate key matching frame_time_map (e.g. 'rep_1', 'rep_2')
        m_rep = re.search(r"rep_?(\d+)", r_dir.name)
        rep_key = f"rep_{m_rep.group(1)}" if m_rep else r_dir.name
        rep_dfs[rep_key] = df_clean

    if not rep_dfs:
        raise FileNotFoundError(
            f"FATAL: Zero valid MM-PBSA CSV deliverables parsed from {mmpbsa_dir}!"
        )

    return rep_dfs


def main():
    args = parse_args()
    system_name = args.system
    temp_k = float(args.temp_k)
    kb = 0.001987204  # Boltzmann constant in kcal/(mol * K)

    labels_p = Path(args.labels).resolve()
    time_map_p = Path(args.time_map).resolve()
    mmpbsa_dir = Path(args.mmpbsa_dir).resolve()

    if not labels_p.exists() or not time_map_p.exists():
        raise FileNotFoundError(
            "FATAL: Required labels array or time map JSON deliverable missing!"
        )

    labels = np.load(labels_p)
    with open(time_map_p, "r", encoding="utf-8") as f:
        time_map = json.load(f)

    if len(labels) != len(time_map):
        raise ValueError(
            f"FATAL: Inconsistent sample sizes! Labels array has {len(labels)}, time map has {len(time_map)}."
        )

    # 1. Index coordinate frames by replicate and time
    df_coords = pd.DataFrame(time_map)
    df_coords["cluster_id"] = labels
    # 1-based sequential frame index within each replicate
    df_coords["rep_frame"] = df_coords.groupby("replicate").cumcount() + 1

    # 2. Load MM-PBSA samples
    rep_energies = load_replicate_mmpbsa(mmpbsa_dir, system_name)

    # 3. Exact 1:1 Bijective Intersection (Zero Heuristics: Failure on Mismatched Grids)
    stratified_samples = []
    for rep_key, df_e in rep_energies.items():
        df_c_rep = df_coords[df_coords["replicate"] == rep_key].copy()
        if df_c_rep.empty:
            raise KeyError(
                f"FATAL SSOT INTEGRITY VIOLATION: Replicate '{rep_key}' present in MM-PBSA "
                "was not found in clustering coordinate map!"
            )

        df_e_rep = df_e.copy()
        # Canonical timestamp key rounded to integration picosecond tolerance
        df_e_rep["join_time"] = (
            df_e_rep["energy_frame"].astype(float)
            if "energy_frame" in df_e_rep
            else np.arange(len(df_e_rep))
        )
        df_c_rep["join_time"] = df_c_rep["time_ns"].round(3)

        n_c = len(df_c_rep)
        n_e = len(df_e_rep)

        # Exact 1:1 length assertion
        if n_c != n_e:
            raise ValueError(
                f"FATAL DATA INTEGRITY VIOLATION for {rep_key}: MM-PBSA evaluated {n_e} samples, "
                f"but clustering evaluated {n_c} frames! Sampling strides have diverged. "
                "Both processes must derive their grid strictly from 'analysis_sampling' in md.yaml."
            )

        df_c_sorted = df_c_rep.sort_values("time_ns").reset_index(drop=True)
        for idx in range(n_e):
            stratified_samples.append(
                {
                    "replicate": rep_key,
                    "time_ns": float(df_c_sorted.at[idx, "time_ns"]),
                    "cluster_id": int(df_c_sorted.at[idx, "cluster_id"]),
                    "delta_g": float(df_e_rep.iloc[idx]["delta_g"]),
                }
            )

    if not stratified_samples:
        raise ValueError(
            "FATAL: Point-to-State projection produced zero valid intersections!"
        )

    df_strat = pd.DataFrame(stratified_samples)
    total_energy_points = len(df_strat)

    # 4. Statistical Moments Calculation per Cluster
    unique_clusters = sorted(list(set(df_strat["cluster_id"])))
    cluster_results = []
    boltzmann_terms = []

    for cid in unique_clusters:
        sub = df_strat[df_strat["cluster_id"] == cid]
        n_points = len(sub)
        pop_fraction = n_points / total_energy_points

        # Replicate-level means to assess between-replicate variance
        rep_means = sub.groupby("replicate")["delta_g"].mean().values
        mean_dg = float(sub["delta_g"].mean())
        std_dg = float(sub["delta_g"].std(ddof=1)) if n_points > 1 else 0.0

        if len(rep_means) > 1:
            sem_dg = float(np.std(rep_means, ddof=1) / np.sqrt(len(rep_means)))
        else:
            sem_dg = float(std_dg / np.sqrt(n_points)) if n_points > 1 else 0.0

        # Term for Boltzmann ensemble aggregate
        boltzmann_weight = pop_fraction * np.exp(-mean_dg / (kb * temp_k))
        boltzmann_terms.append(boltzmann_weight)

        cluster_results.append(
            {
                "cluster_id": int(cid),
                "sample_count": n_points,
                "population_fraction": round(pop_fraction, 4),
                "delta_g_mean": round(mean_dg, 3),
                "delta_g_sem": round(sem_dg, 3),
                "delta_g_std": round(std_dg, 3),
            }
        )

    # 5. Master Boltzmann Ensemble Free Energy
    sum_weights = sum(boltzmann_terms)
    if sum_weights > 0:
        boltzmann_total_dg = float(-kb * temp_k * np.log(sum_weights))
    else:
        boltzmann_total_dg = float(df_strat["delta_g"].mean())

    summary_manifest = {
        "system": system_name,
        "temperature_kelvin": temp_k,
        "total_independent_samples": total_energy_points,
        "boltzmann_total_delta_g": round(boltzmann_total_dg, 3),
        "clusters": cluster_results,
    }

    with open(args.out_json, "w", encoding="utf-8") as jf:
        json.dump(summary_manifest, jf, indent=2)

    # 6. Render Publication Multi-Replica Color-Coded Trajectory Plot
    palette = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
    ]
    color_map = {
        cid: palette[idx % len(palette)] for idx, cid in enumerate(unique_clusters)
    }

    unique_reps = sorted(list(set(df_coords["replicate"])))
    n_reps = len(unique_reps)

    n_cols = 2 if n_reps > 1 else 1
    n_rows = int(np.ceil(n_reps / float(n_cols)))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(14, max(4, n_rows * 2.5)),
        dpi=250,
        sharex=True,
        sharey=True,
    )
    if n_reps == 1:
        axes_flat = [axes]
    else:
        axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

    # Look for RMSD tables from md_metrics
    metrics_base = (
        Path(args.metrics_dir).resolve()
        if args.metrics_dir
        else Path.cwd() / "md_runs" / system_name
    )

    for idx, r_key in enumerate(unique_reps):
        ax = axes_flat[idx]
        df_c_r = df_coords[df_coords["replicate"] == r_key].sort_values("time_ns")
        times = df_c_r["time_ns"].values
        c_ids = df_c_r["cluster_id"].values

        # Locate clean RMSD table for this replicate if available
        m_num = re.search(r"(\d+)", r_key)
        rep_num = m_num.group(1) if m_num else "1"
        rmsd_file = metrics_base / f"rep_{rep_num}" / "metrics" / "rmsd.csv"

        has_rmsd = False
        if rmsd_file.exists():
            try:
                df_rmsd = pd.read_csv(rmsd_file)
                lig_col = next((c for c in df_rmsd.columns if "LIG" in c.upper()), None)
                if lig_col and len(df_rmsd) >= len(times) // 2:
                    y_vals = np.interp(times, df_rmsd["Time_ns"], df_rmsd[lig_col])
                    y_label = "Ligand RMSD (Å)"
                    has_rmsd = True
            except Exception:
                pass

        if not has_rmsd:
            # Fallback: plot instantaneous Delta G values for this replicate
            sub_e = df_strat[df_strat["replicate"] == r_key]
            times = sub_e["time_ns"].values
            y_vals = sub_e["delta_g"].values
            c_ids = sub_e["cluster_id"].values
            y_label = "ΔG (kcal/mol)"

        if len(times) > 1 and len(y_vals) > 1:
            points = np.array([times, y_vals]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            seg_colors = [color_map.get(cid, "gray") for cid in c_ids[:-1]]

            lc = LineCollection(segments, colors=seg_colors, linewidth=1.8)
            ax.add_collection(lc)
            ax.set_xlim(times.min(), times.max())
            ax.set_ylim(
                min(0.0, float(np.min(y_vals)) * 0.9), float(np.max(y_vals)) * 1.1
            )

        ax.set_title(f"{r_key}", fontsize=10, fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.5)

    # Hide unused subplot panels
    for j in range(n_reps, len(axes_flat)):
        axes_flat[j].axis("off")

    # Legend reflecting thermodynamics
    custom_lines = []
    legend_labels = []
    for c_stat in cluster_results:
        cid = c_stat["cluster_id"]
        c_col = color_map[cid]
        line = plt.Line2D([0], [0], color=c_col, lw=2.5)
        lbl = (
            f"Basin #{cid}: ΔG = {c_stat['delta_g_mean']:.2f} ± {c_stat['delta_g_sem']:.2f} kcal/mol "
            f"({c_stat['population_fraction']*100:.1f}%)"
        )
        custom_lines.append(line)
        legend_labels.append(lbl)

    fig.legend(
        custom_lines,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=min(3, len(cluster_results)),
        frameon=True,
        fontsize=10,
    )

    fig.text(
        0.5, 0.01, "Simulation Time (ns)", ha="center", fontsize=11, fontweight="bold"
    )
    fig.text(
        0.02,
        0.5,
        y_label,
        va="center",
        rotation="vertical",
        fontsize=11,
        fontweight="bold",
    )
    plt.tight_layout(rect=[0.03, 0.03, 0.98, 0.95])
    plt.savefig(args.out_plot)
    plt.close()

    print(
        f"✓ Stratification completed successfully -> {args.out_json} & {args.out_plot}",
        flush=True,
    )


if __name__ == "__main__":
    main()
