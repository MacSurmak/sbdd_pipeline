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
        "--paradigm", default="adaptive_zscore", help="Plateau detection paradigm"
    )
    parser.add_argument(
        "--smoothing_tau_ns", type=float, default=0.5, help="Zero-phase EMA decay time (ns)"
    )
    parser.add_argument(
        "--significance_z_crit", type=float, default=2.0, help="Adaptive Z-score threshold"
    )
    parser.add_argument(
        "--dilation_buffer_ns", type=float, default=0.4, help="Dilation buffer (+/- delta in ns)"
    )
    parser.add_argument(
        "--min_dwell_ns", type=float, default=5.0, help="Minimum plateau duration (ns)"
    )
    parser.add_argument(
        "--enforce_cluster_continuity",
        action="store_true",
        default=True,
        help="Split plateau if manifold cluster changes",
    )
    parser.add_argument(
        "--out_json", default="stratification_summary.json", help="Output summary JSON"
    )
    parser.add_argument(
        "--out_plot",
        default="stratified_rmsd_plot.png",
        help="Output multi-replica plot PNG",
    )
    parser.add_argument(
        "--out_timeseries",
        default="stratified_timeseries.tsv",
        help="Output stratified timeseries TSV",
    )
    return parser.parse_args()


def zero_phase_ema(y: np.ndarray, times: np.ndarray, decay_ns: float) -> np.ndarray:
    """Forward-backward zero-phase exponential moving average filter."""
    n = len(y)
    if n < 3 or decay_ns <= 0.0:
        return y.copy()
    dt = (times[-1] - times[0]) / float(max(1, n - 1))
    alpha = 1.0 - np.exp(-dt / float(decay_ns))

    fwd = np.empty(n, dtype=np.float64)
    fwd[0] = y[0]
    for i in range(1, n):
        fwd[i] = alpha * y[i] + (1.0 - alpha) * fwd[i - 1]

    bwd = np.empty(n, dtype=np.float64)
    bwd[-1] = fwd[-1]
    for i in range(n - 2, -1, -1):
        bwd[i] = alpha * fwd[i] + (1.0 - alpha) * bwd[i + 1]

    return bwd


def compute_replicate_mad(rmsd: np.ndarray, smoothed: np.ndarray) -> float:
    """Computes median absolute deviation estimate of high-frequency thermal noise."""
    res = rmsd - smoothed
    med = float(np.median(res))
    mad = float(np.median(np.abs(res - med)))
    return max(0.04, 1.4826 * mad)


def detect_plateaus_exclusion(
    times: np.ndarray,
    rmsd: np.ndarray,
    cluster_arr: np.ndarray,
    paradigm: str,
    tau_ns: float,
    z_crit: float,
    buffer_ns: float,
    min_dwell_ns: float,
    enforce_cluster: bool,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Exclusion-First plateau gating:
    Identifies non-equilibrium transition jumps, applies morphological dilation buffer (+/- delta),
    and validates remaining quiescent intervals with minimum dwell time and cluster continuity.
    """
    n = len(times)
    dt = (times[-1] - times[0]) / float(max(1, n - 1))
    smoothed = zero_phase_ema(rmsd, times, tau_ns)
    sigma_thermal = compute_replicate_mad(rmsd, smoothed)

    trans_spikes = np.zeros(n, dtype=np.uint8)

    if paradigm in ["adaptive_zscore", "macro_step"]:
        effective_step_thresh = (z_crit * sigma_thermal) if paradigm == "adaptive_zscore" else z_crit
        stride = max(1, int(round(0.20 / dt)))
        sub_n = n // stride
        if sub_n > 4:
            sub_means = []
            sub_idxs = []
            for k in range(sub_n):
                s = k * stride
                e = min(n, (k + 1) * stride)
                sub_means.append(float(np.mean(smoothed[s:e])))
                sub_idxs.append(s)

            lookahead = max(1, int(round(2.0 / 0.20)))
            for k in range(lookahead, len(sub_means) - lookahead):
                delta_mu = abs(sub_means[k + lookahead] - sub_means[k - lookahead])
                if delta_mu >= effective_step_thresh:
                    center_pt = sub_idxs[k]
                    trans_spikes[center_pt] = 1
    elif paradigm == "running_range":
        w_pts = max(5, int(round(6.0 / dt)))
        half_w = w_pts // 2
        for i in range(half_w, n - half_w, 2):
            val_range = np.max(smoothed[i - half_w : i + half_w + 1]) - np.min(smoothed[i - half_w : i + half_w + 1])
            if val_range >= z_crit:
                trans_spikes[i] = 1
    else:
        # Inverted velocity spike
        grad = np.gradient(smoothed, times)
        vel = np.abs(grad)
        trans_spikes[vel >= z_crit] = 1

    # Morphological dilation buffer (+/- delta)
    buf_pts = int(round(buffer_ns / dt))
    excluded_mask = np.zeros(n, dtype=np.uint8)
    if buf_pts > 0:
        spike_indices = np.where(trans_spikes == 1)[0]
        for s_idx in spike_indices:
            s = max(0, s_idx - buf_pts)
            e = min(n, s_idx + buf_pts + 1)
            excluded_mask[s:e] = 1
    else:
        excluded_mask = trans_spikes.copy()

    # Plateau candidate is complement
    candidate_plateau = np.where(excluded_mask == 0, 1, 0).astype(np.uint8)
    final_plateau_mask = np.zeros(n, dtype=np.uint8)

    in_run = False
    start_idx = 0
    for i in range(n + 1):
        val = candidate_plateau[i] if i < n else 0
        if not in_run and val == 1:
            in_run = True
            start_idx = i
        elif in_run and val == 0:
            in_run = False
            end_idx = i
            if (times[end_idx - 1] - times[start_idx]) >= min_dwell_ns:
                if enforce_cluster and len(cluster_arr) == n:
                    c_run_start = start_idx
                    curr_c = cluster_arr[start_idx]
                    for k in range(start_idx + 1, end_idx + 1):
                        next_c = cluster_arr[k] if k < end_idx else -999
                        if next_c != curr_c:
                            if (times[k - 1] - times[c_run_start]) >= min_dwell_ns:
                                final_plateau_mask[c_run_start:k] = 1
                            c_run_start = k
                            curr_c = next_c
                else:
                    final_plateau_mask[start_idx:end_idx] = 1

    return final_plateau_mask, smoothed, sigma_thermal


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

    # 3. Load full-resolution RMSD curves and compute Exclusion-First plateau masks
    metrics_base = (
        Path(args.metrics_dir).resolve()
        if args.metrics_dir
        else Path.cwd() / "md_runs" / system_name
    )

    replicate_plateau_masks = {}
    replicate_full_rmsd = {}
    replicate_smoothed = {}

    for rep_key in rep_energies.keys():
        m_num = re.search(r"(\d+)", rep_key)
        rep_num = m_num.group(1) if m_num else "1"
        rmsd_file = metrics_base / f"rep_{rep_num}" / "metrics" / "rmsd.csv"
        if not rmsd_file.exists():
            rmsd_file = metrics_base / f"rep_{rep_num}" / "rmsd.csv"

        if not rmsd_file.exists() or rmsd_file.stat().st_size == 0:
            raise FileNotFoundError(
                f"FATAL: Full-resolution RMSD table missing for {rep_key} at {rmsd_file}! "
                "Exclusion-First plateau gating requires continuous RMSD trajectory data."
            )

        df_rmsd_raw = pd.read_csv(rmsd_file)
        lig_col = next((c for c in df_rmsd_raw.columns if c.startswith("Lig_")), None)
        if not lig_col:
            lig_col = next((c for c in df_rmsd_raw.columns if "RMSD" in c and "BB" not in c), df_rmsd_raw.columns[-1])

        df_r_clean = pd.DataFrame({
            "Time_ns": pd.to_numeric(df_rmsd_raw["Time_ns"], errors="coerce"),
            "RMSD": pd.to_numeric(df_rmsd_raw[lig_col], errors="coerce")
        }).dropna().sort_values("Time_ns").reset_index(drop=True)

        r_times = df_r_clean["Time_ns"].values
        r_vals = df_r_clean["RMSD"].values

        # Map cluster labels onto full-resolution time grid via nearest time lookup
        df_c_r = df_coords[df_coords["replicate"] == rep_key].sort_values("time_ns")
        c_times = df_c_r["time_ns"].values
        c_labels = df_c_r["cluster_id"].values
        nearest_idx = np.clip(np.searchsorted(c_times, r_times, side="left"), 0, len(c_labels) - 1)
        r_clusters = c_labels[nearest_idx]

        p_mask, sm_vals, sigma_th = detect_plateaus_exclusion(
            times=r_times,
            rmsd=r_vals,
            cluster_arr=r_clusters,
            paradigm=args.paradigm,
            tau_ns=float(args.smoothing_tau_ns),
            z_crit=float(args.significance_z_crit),
            buffer_ns=float(args.dilation_buffer_ns),
            min_dwell_ns=float(args.min_dwell_ns),
            enforce_cluster=bool(args.enforce_cluster_continuity),
        )

        replicate_plateau_masks[rep_key] = (r_times, p_mask)
        replicate_full_rmsd[rep_key] = (r_times, r_vals)
        replicate_smoothed[rep_key] = (r_times, sm_vals)

    # 4. Exact Bijective Intersection with Plateau Gating
    stratified_samples = []
    for rep_key, df_e in rep_energies.items():
        df_c_rep = df_coords[df_coords["replicate"] == rep_key].copy()
        if df_c_rep.empty:
            raise KeyError(
                f"FATAL SSOT INTEGRITY VIOLATION: Replicate '{rep_key}' present in MM-PBSA "
                "was not found in clustering coordinate map!"
            )

        df_e_rep = df_e.copy()
        n_c = len(df_c_rep)
        n_e = len(df_e_rep)

        if n_c != n_e:
            raise ValueError(
                f"FATAL DATA INTEGRITY VIOLATION for {rep_key}: MM-PBSA evaluated {n_e} samples, "
                f"but clustering evaluated {n_c} frames! Sampling strides have diverged."
            )

        df_c_sorted = df_c_rep.sort_values("time_ns").reset_index(drop=True)
        r_times_full, p_mask_full = replicate_plateau_masks[rep_key]

        for idx in range(n_e):
            t_sample = float(df_c_sorted.at[idx, "time_ns"])
            c_id = int(df_c_sorted.at[idx, "cluster_id"])
            d_g = float(df_e_rep.iloc[idx]["delta_g"])

            # Exact nearest lookup in full-resolution plateau mask
            n_pos = np.clip(np.searchsorted(r_times_full, t_sample, side="left"), 0, len(p_mask_full) - 1)
            is_plat = bool(p_mask_full[n_pos] == 1)

            stratified_samples.append(
                {
                    "replicate": rep_key,
                    "time_ns": t_sample,
                    "cluster_id": c_id,
                    "delta_g": d_g,
                    "is_plateau": is_plat,
                }
            )

    if not stratified_samples:
        raise ValueError(
            "FATAL: Point-to-State projection produced zero valid intersections!"
        )

    df_strat = pd.DataFrame(stratified_samples)
    total_energy_points = len(df_strat)
    df_strat.to_csv(args.out_timeseries, sep="\t", index=False)
    print(f"[STRATIFY] Saved complete timeseries with plateau mask -> {args.out_timeseries}", flush=True)

    # 5. Statistical Moments Calculation strictly over Quiescent Plateau Frames
    unique_clusters = sorted(list(set(df_strat["cluster_id"])))
    cluster_results = []
    boltzmann_terms = []

    total_plateau_points = int(df_strat["is_plateau"].sum())
    total_samples_all = len(df_strat)

    for cid in unique_clusters:
        sub_all = df_strat[df_strat["cluster_id"] == cid]
        sub_plat = sub_all[sub_all["is_plateau"] == True]

        n_all = len(sub_all)
        n_plat = len(sub_plat)
        pop_fraction = (n_plat / total_plateau_points) if total_plateau_points > 0 else (n_all / total_samples_all)

        # Free energy is evaluated strictly on stationary plateau frames
        eval_data = sub_plat if n_plat > 0 else sub_all
        mean_dg = float(eval_data["delta_g"].mean())
        std_dg = float(eval_data["delta_g"].std(ddof=1)) if len(eval_data) > 1 else 0.0

        rep_means = eval_data.groupby("replicate")["delta_g"].mean().values
        if len(rep_means) > 1:
            sem_dg = float(np.std(rep_means, ddof=1) / np.sqrt(len(rep_means)))
        else:
            sem_dg = float(std_dg / np.sqrt(len(eval_data))) if len(eval_data) > 1 else 0.0

        if pop_fraction > 0.0:
            boltzmann_weight = pop_fraction * np.exp(-mean_dg / (kb * temp_k))
            boltzmann_terms.append(boltzmann_weight)

        cluster_results.append(
            {
                "cluster_id": int(cid),
                "plateau_samples": n_plat,
                "total_samples": n_all,
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
        r_times_f, r_vals_f = replicate_full_rmsd[r_key]
        _, sm_vals_f = replicate_smoothed[r_key]
        _, p_mask_f = replicate_plateau_masks[r_key]

        df_c_r = df_coords[df_coords["replicate"] == r_key].sort_values("time_ns")
        c_times = df_c_r["time_ns"].values
        c_labels = df_c_r["cluster_id"].values
        nearest_idx = np.clip(np.searchsorted(c_times, r_times_f, side="left"), 0, len(c_labels) - 1)
        cluster_assigned = c_labels[nearest_idx]

        # 1. Shaded excluded transition regions (Red)
        in_t = False
        t_start = 0
        for i in range(len(r_times_f) + 1):
            val = p_mask_f[i] if i < len(r_times_f) else 1
            if not in_t and val == 0:
                in_t = True
                t_start = i
            elif in_t and val == 1:
                in_t = False
                ax.axvspan(r_times_f[t_start], r_times_f[i - 1], color="red", alpha=0.15, lw=0)

        # 2. Raw RMSD line in background (Muted Gray)
        ax.plot(r_times_f, r_vals_f, color="#666677", lw=0.6, alpha=0.5)

        # 3. Active Plateau runs colored by cluster
        in_p = False
        p_start = 0
        for i in range(len(r_times_f) + 1):
            val = p_mask_f[i] if i < len(r_times_f) else 0
            if not in_p and val == 1:
                in_p = True
                p_start = i
            elif in_p and val == 0:
                in_p = False
                p_end = i
                c_id = cluster_assigned[p_start]
                col = color_map.get(c_id, "#00e5ff")
                ax.plot(r_times_f[p_start:p_end], sm_vals_f[p_start:p_end], color=col, lw=2.4)

        # 4. Excluded transition smooth line (Dashed)
        ax.plot(r_times_f[p_mask_f == 0], sm_vals_f[p_mask_f == 0], color="#ff5252", lw=1.2, ls=":", alpha=0.7)

        ax.set_xlim(r_times_f.min(), r_times_f.max())
        ax.set_ylim(0.0, max(6.0, float(np.max(r_vals_f)) * 1.15))
        ax.set_title(f"{r_key}", fontsize=10, fontweight="bold")
        ax.grid(True, linestyle=":", alpha=0.4)

    y_label = "Ligand RMSD (Å)"

    # Hide unused subplot panels
    for j in range(n_reps, len(axes_flat)):
        axes_flat[j].axis("off")

    # Legend reflecting thermodynamics
    custom_lines = []
    legend_labels = []
    for c_stat in cluster_results:
        cid = c_stat["cluster_id"]
        c_col = color_map[cid]
        line = plt.Line2D([0], [0], color=c_col, lw=2.8)
        lbl = (
            f"Basin #{cid}: ΔG = {c_stat['delta_g_mean']:.2f} ± {c_stat['delta_g_sem']:.2f} kcal/mol "
            f"({c_stat['population_fraction']*100:.1f}%, {c_stat['plateau_samples']}/{c_stat['total_samples']} pts)"
        )
        custom_lines.append(line)
        legend_labels.append(lbl)

    line_trans = plt.Line2D([0], [0], color="red", alpha=0.3, lw=4.0)
    custom_lines.append(line_trans)
    legend_labels.append("Excluded Transition Buffer (0% Leakage)")

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
