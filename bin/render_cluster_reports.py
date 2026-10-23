#!/usr/bin/env python3
"""
bin/render_cluster_reports.py
Atomic Stage 4: Boltzmann Free Energy Surface (FES) calculation,
Plotly HTML interactive dashboard, publication PNG, and portable relative PyMOL session.
"""

import argparse
import json
import shutil
import sys
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def compute_fes_2d(c1: np.ndarray, c2: np.ndarray, temp_k: float, bins: int = 80):
    """
    Computes physical Delta G(x, y) = -kB * T * ln(P(x, y) / Pmax) in kcal/mol.
    Strictly masks unvisited phase space as NaN to render empty background
    without generating artificial energy plateaus.
    """
    from scipy.ndimage import gaussian_filter

    kb = 0.001987204  # Boltzmann constant in kcal / (mol * K)
    hist, xedges, yedges = np.histogram2d(c1, c2, bins=bins, density=True)
    hist = hist.T

    smoothed_hist = gaussian_filter(hist, sigma=1.0)
    max_val = np.max(smoothed_hist)
    if max_val <= 0:
        raise ValueError("FATAL: Maximum probability density for FES is zero!")

    prob = smoothed_hist / max_val
    # Strictly mask empty phase space (prob < 1e-4) with NaN
    empty_mask = prob < 1e-4
    prob_clean = np.where(empty_mask, np.nan, prob)

    with np.errstate(invalid="ignore"):
        dg = -kb * temp_k * np.log(prob_clean)
        dg -= np.nanmin(dg)

    xc = 0.5 * (xedges[:-1] + xedges[1:])
    yc = 0.5 * (yedges[:-1] + yedges[1:])
    return xc, yc, dg


def main():
    parser = argparse.ArgumentParser(description="Render thermodynamic reports and portable PyMOL sessions")
    parser.add_argument("--struct", required=True, help="Clean reference structure PDB")
    parser.add_argument("--umap", required=True, help="Input umap_2d.npy array")
    parser.add_argument("--projections", required=True, help="Input projections.npy array")
    parser.add_argument("--labels", required=True, help="Input cluster_labels.npy array")
    parser.add_argument("--medoids_json", required=True, help="Input medoids_summary.json deliverable")
    parser.add_argument("--time_map", required=True, help="Input frame_time_map.json mapping")
    parser.add_argument("--feature_meta", required=True, help="Input feature_meta.json deliverable")
    parser.add_argument("--align_group", required=True, help="Alignment basis group name")
    parser.add_argument("--cluster_group", required=True, help="Clustering feature group name")
    parser.add_argument("--dim_reduction", required=True, help="Dimensionality reduction method (pca or tica)")
    parser.add_argument("--system", required=True, help="Canonical parent system name")
    parser.add_argument("--extracted_json", required=True, help="Input extracted_manifest.json with structure deliverables")
    parser.add_argument("--temp_k", type=float, required=True, help="Simulation temperature in Kelvin")
    parser.add_argument("--outdir", required=True, help="Output directory for reports")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    with open(Path(args.feature_meta).resolve(), "r", encoding="utf-8") as f:
        feat_meta = json.load(f)

    umap_2d = np.load(Path(args.umap).resolve())
    proj = np.load(Path(args.projections).resolve())
    labels = np.load(Path(args.labels).resolve())

    with open(Path(args.medoids_json).resolve(), "r", encoding="utf-8") as f:
        medoids = json.load(f)
    with open(Path(args.time_map).resolve(), "r", encoding="utf-8") as f:
        time_map = json.load(f)

    n_frames = umap_2d.shape[0]
    palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"]

    # Copy reference PDB into outdir for 100% self-contained relative PyMOL sessions
    ref_pdb = outdir / "reference.pdb"
    shutil.copy2(Path(args.struct).resolve(), ref_pdb)

    # 1. Compute Free Energy Surfaces
    u_xc, u_yc, u_dg = compute_fes_2d(umap_2d[:, 0], umap_2d[:, 1], args.temp_k)
    p_xc, p_yc, p_dg = compute_fes_2d(proj[:, 0], proj[:, 1], args.temp_k)

    # 2. Static Publication Figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=250)

    ax1 = axes[0, 0]
    cp1 = ax1.contourf(u_xc, u_yc, u_dg, levels=18, cmap="viridis")
    cbar1 = fig.colorbar(cp1, ax=ax1)
    cbar1.set_label("ΔG (kcal/mol)")
    for m in medoids:
        f = m["medoid_frame"]
        ax1.scatter(umap_2d[f, 0], umap_2d[f, 1], marker="*", s=180, color="red", edgecolors="black", zorder=5)
        ax1.text(umap_2d[f, 0] + 0.2, umap_2d[f, 1] + 0.2, f"#{m['rank']} ({m['population_pct']}%)", color="white", fontweight="bold", fontsize=9)
    ax1.set_title("1. Free Energy Surface ΔG on UMAP (kcal/mol)", fontweight="bold")
    ax1.set_xlabel("UMAP 1")
    ax1.set_ylabel("UMAP 2")

    ax2 = axes[0, 1]
    for m in medoids:
        cid = m["raw_cluster_id"]
        mask = labels == cid
        col = palette[(m["rank"] - 1) % len(palette)]
        ax2.scatter(umap_2d[mask, 0], umap_2d[mask, 1], s=3, alpha=0.5, color=col, label=f"Basin #{m['rank']} ({m['population_pct']}%)")
    ax2.set_title(f"2. Partitioned Macro-Basins (K={len(medoids)})", fontweight="bold")
    ax2.set_xlabel("UMAP 1")
    ax2.set_ylabel("UMAP 2")
    ax2.legend(markerscale=4, frameon=True)

    ax3 = axes[1, 0]
    cp3 = ax3.contourf(p_xc, p_yc, p_dg, levels=18, cmap="viridis")
    cbar3 = fig.colorbar(cp3, ax=ax3)
    cbar3.set_label("ΔG (kcal/mol)")
    for m in medoids:
        f = m["medoid_frame"]
        ax3.scatter(proj[f, 0], proj[f, 1], marker="*", s=180, color="red", edgecolors="black", zorder=5)
        ax3.text(proj[f, 0] + 0.8, proj[f, 1] + 0.8, f"#{m['rank']}", color="white", fontweight="bold", fontsize=9)
    ax3.set_title("3. Free Energy Surface ΔG on Component Space", fontweight="bold")
    ax3.set_xlabel("Component 1")
    ax3.set_ylabel("Component 2")

    ax4 = axes[1, 1]
    unique_reps = sorted(list(set(item["replicate"] for item in time_map)))
    bar_width = 0.5
    bottom_vals = np.zeros(len(medoids))
    for r in unique_reps:
        counts = [m["replicate_breakdown"].get(r, 0) for m in medoids]
        ax4.bar([f"Basin #{m['rank']}" for m in medoids], counts, bottom=bottom_vals, label=r, width=bar_width)
        bottom_vals += np.array(counts)
    ax4.set_title("4. Replicate Breakdown per Basin", fontweight="bold")
    ax4.set_ylabel("Frame Count")
    ax4.legend(frameon=True)

    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    png_path = outdir / "macro_clustering_landscape.png"
    plt.savefig(png_path)
    plt.close()

    # 3. Interactive Plotly Dashboard
    step = max(1, n_frames // 10000)
    sub = slice(0, n_frames, step)

    fig_html = make_subplots(rows=1, cols=2, subplot_titles=("UMAP Basins & Energy Medoids", "Component Representation of Basins"))
    for m in medoids:
        cid = m["raw_cluster_id"]
        mask = labels == cid
        sub_mask = mask[sub]
        col = palette[(m["rank"] - 1) % len(palette)]
        fig_html.add_trace(
            go.Scatter(x=umap_2d[sub][sub_mask, 0], y=umap_2d[sub][sub_mask, 1], mode="markers",
                       name=f"Basin #{m['rank']} ({m['population_pct']}%)", marker=dict(size=3, color=col, opacity=0.6)),
            row=1, col=1
        )
        fig_html.add_trace(
            go.Scatter(x=proj[sub][sub_mask, 0], y=proj[sub][sub_mask, 1], mode="markers",
                       name=f"Basin #{m['rank']} ({m['population_pct']}%)", marker=dict(size=3, color=col, opacity=0.6), showlegend=False),
            row=1, col=2
        )

    for m in medoids:
        f = m["medoid_frame"]
        fig_html.add_trace(
            go.Scatter(x=[umap_2d[f, 0]], y=[umap_2d[f, 1]], mode="markers+text",
                       text=[f"#{m['rank']}"], textposition="top center",
                       marker=dict(size=14, color="red", symbol="star", line=dict(color="black", width=1.5)),
                       name=f"Medoid #{m['rank']}"),
            row=1, col=1
        )
        fig_html.add_trace(
            go.Scatter(x=[proj[f, 0]], y=[proj[f, 1]], mode="markers+text",
                       text=[f"#{m['rank']}"], textposition="top center",
                       marker=dict(size=14, color="red", symbol="star", line=dict(color="black", width=1.5)),
                       showlegend=False),
            row=1, col=2
        )

    fig_html.update_layout(height=650, width=1300, template="plotly_dark", title_text="Conformational Macro-Basin Landscape Dashboard")
    html_path = outdir / "macro_clustering_interactive.html"
    fig_html.write_html(str(html_path))

    # 4. Portable PyMOL Session Script (STRICTLY RELATIVE PATHS)
    pml_lines = [
        "# PyMOL Centroid Inspection Session (Fully Portable)",
        "reinitialize",
        "load reference.pdb, reference_pdb",
        "color gray70, reference_pdb",
        "show cartoon, reference_pdb",
        "set cartoon_transparency, 0.4, reference_pdb",
        "",
    ]
    pymol_colors = ["cyan", "orange", "palegreen", "salmon", "magenta", "yellow", "lightblue", "lightpink", "wheat"]
    for m in medoids:
        rank = m["rank"]
        pop = m["population_pct"]
        col = pymol_colors[(rank - 1) % len(pymol_colors)]
        obj_name = f"centroid_{rank}_{pop:.1f}pct"

        pml_lines.extend([
            f"load centroid_{rank}_complex.pdb, {obj_name}",
            f"align {obj_name} and name CA, reference_pdb and name CA",
            f"color {col}, {obj_name}",
            f"show cartoon, {obj_name}",
            f"show sticks, {obj_name} and not protein and not (resname HOH or resname TIP3 or resname WAT)",
            "",
        ])
    pml_lines.extend(["zoom reference_pdb, buffer=5.0", "util.cnc"])

    pml_path = outdir / "inspect_clusters.pml"
    pml_path.write_text("\n".join(pml_lines) + "\n")

    with open(Path(args.extracted_json).resolve(), "r", encoding="utf-8") as ef:
        extracted_structures = json.load(ef)

    manifest = {
        "system": args.system,
        "n_frames": n_frames,
        "n_clusters": len(medoids),
        "align_group": args.align_group,
        "cluster_group": args.cluster_group,
        "dim_reduction": args.dim_reduction,
        "html_report": html_path.name,
        "png_overview": png_path.name,
        "pymol_session": pml_path.name,
        "basins": medoids,
        "structures": extracted_structures,
        "feature_metadata": feat_meta,
    }
    with open(outdir / "clusters_manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"✓ Generated reports and portable PyMOL session in {outdir}", flush=True)


if __name__ == "__main__":
    main()
