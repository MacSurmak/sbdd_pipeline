#!/usr/bin/env python3
"""
bin/summarize_results.py
Master aggregator:
Parses Glide CSVs, GNINA Dock SDFs, GNINA Rescore SDFs, Prime MM-GBSA CSVs,
and outputs complete, structured analysis CSV tables:
1. ranking.csv (Consensus ranked table with Ultimate_Score, Percentile Ranks, Glide, MMGBSA, CNN Affinity, Pose RMSD)
2. cross_frame_matrix.csv (Pivot matrix of Ultimate_Score across all frames)
3. summary.json (JSON dump of ranking)
"""
import argparse
import json
from pathlib import Path
import sys
from typing import Final
import numpy as np
import pandas as pd
from rdkit import Chem

SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
FRAMEWORK_DIR: Final[Path] = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(FRAMEWORK_DIR) not in sys.path:
    sys.path.insert(0, str(FRAMEWORK_DIR))

from cli.core.chem import compute_symmetry_rmsd
from utils.parsers import parse_glide_csv, parse_gnina_sdf, parse_prime_csv


def compute_pose_consensus_rmsd(df: pd.DataFrame, dock_dir: Path) -> pd.DataFrame:
    """Calculates in-situ symmetry-aware heavy-atom RMSD between top Glide and Uni-Dock poses."""

    def _index_top_poses(filepath: Path) -> dict[str, Chem.Mol]:
        poses_by_title: dict[str, Chem.Mol] = {}
        if not filepath.exists() or filepath.stat().st_size == 0:
            return poses_by_title
        suppl = Chem.SDMolSupplier(str(filepath), removeHs=False, sanitize=False)
        for m in suppl:
            if m is not None:
                m_title = m.GetProp("_Name").strip() if m.HasProp("_Name") else ""
                if not m_title and m.HasProp("s_m_title"):
                    m_title = m.GetProp("s_m_title").strip()
                if m_title and m_title not in poses_by_title:
                    poses_by_title[m_title] = m
        return poses_by_title

    glide_cache: dict[str, dict[str, Chem.Mol]] = {}
    unidock_cache: dict[str, dict[str, Chem.Mol]] = {}
    for rec_val in df["Receptor"].dropna().unique():
        r_id = str(rec_val).strip()
        glide_sdf = (
            dock_dir / f"{r_id}_glide_rescored.sdf"
            if (dock_dir / f"{r_id}_glide_rescored.sdf").exists()
            else dock_dir / f"{r_id}_dock.sdf"
        )
        unidock_sdf = (
            dock_dir / f"{r_id}_unidock_rescored.sdf"
            if (dock_dir / f"{r_id}_unidock_rescored.sdf").exists()
            else dock_dir / f"{r_id}_unidock_dock.sdf"
        )
        glide_cache[r_id] = _index_top_poses(glide_sdf)
        unidock_cache[r_id] = _index_top_poses(unidock_sdf)

    rmsd_list = []
    for _, row in df.iterrows():
        rec_id = str(row["Receptor"]).strip()
        title = str(row["Title"]).strip()

        mg = (
            glide_cache[rec_id][title]
            if (rec_id in glide_cache and title in glide_cache[rec_id])
            else None
        )
        mu = (
            unidock_cache[rec_id][title]
            if (rec_id in unidock_cache and title in unidock_cache[rec_id])
            else None
        )
        if mg is not None and mu is not None:
            val = compute_symmetry_rmsd(mg, mu)
            rmsd_list.append(round(val, 3))
        else:
            rmsd_list.append(None)

    df["Pose_RMSD"] = rmsd_list
    return df


def calculate_consensus_score(df: pd.DataFrame) -> pd.DataFrame:
    """
    Computes 3-axis Orthogonal Percentile Consensus Score (0-100):
      1. PR_Physics: mean percentile rank of native physical energies (GlideScore, UniDock_Vina, MMGBSA).
      2. PR_Validity: mean percentile rank of 3D-CNN pose plausibility (Glide_CNNscore, UniDock_CNNscore).
      3. PR_Pose_RMSD: percentile rank of cross-engine geometric convergence (Pose_RMSD, lower is better).
    """
    phys_cols: list[str] = []
    if "GlideScore" in df.columns and df["GlideScore"].notna().any():
        df["PR_Glide"] = df["GlideScore"].rank(ascending=False, pct=True).fillna(0.0)
        phys_cols.append("PR_Glide")
    if "UniDock_Vina" in df.columns and df["UniDock_Vina"].notna().any():
        df["PR_UniDock"] = df["UniDock_Vina"].rank(ascending=False, pct=True).fillna(0.0)
        phys_cols.append("PR_UniDock")
    if "MMGBSA" in df.columns and df["MMGBSA"].notna().any():
        df["PR_MMGBSA"] = df["MMGBSA"].rank(ascending=False, pct=True).fillna(0.0)
        phys_cols.append("PR_MMGBSA")

    val_cols: list[str] = []
    if "Glide_CNNscore" in df.columns and df["Glide_CNNscore"].notna().any():
        df["PR_Glide_CNNscore"] = (
            df["Glide_CNNscore"].rank(ascending=True, pct=True).fillna(0.0)
        )
        val_cols.append("PR_Glide_CNNscore")
    if "UniDock_CNNscore" in df.columns and df["UniDock_CNNscore"].notna().any():
        df["PR_UniDock_CNNscore"] = (
            df["UniDock_CNNscore"].rank(ascending=True, pct=True).fillna(0.0)
        )
        val_cols.append("PR_UniDock_CNNscore")

    axis_cols: list[str] = []
    if phys_cols:
        df["PR_Physics"] = df[phys_cols].sum(axis=1) / float(len(phys_cols))
        axis_cols.append("PR_Physics")
    if val_cols:
        df["PR_Validity"] = df[val_cols].sum(axis=1) / float(len(val_cols))
        axis_cols.append("PR_Validity")
    if "Pose_RMSD" in df.columns and df["Pose_RMSD"].notna().any():
        df["PR_Pose_RMSD"] = df["Pose_RMSD"].rank(ascending=False, pct=True).fillna(0.0)
        axis_cols.append("PR_Pose_RMSD")

    if not axis_cols:
        df["Ultimate_Score"] = 0.0
    else:
        df["Ultimate_Score"] = (df[axis_cols].sum(axis=1) / float(len(axis_cols))) * 100.0
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate project results into master ranking"
    )
    parser.add_argument(
        "--docking_dir",
        required=True,
        help="Directory containing docking poses and tables",
    )
    parser.add_argument(
        "--plan",
        default=None,
        help="Explicit JSON execution plan (dock_plan.json or rescore_plan.json)",
    )
    parser.add_argument("--analysis_dir", help="Directory containing MM-GBSA csv files")
    parser.add_argument(
        "--native_ligand", default="", help="Native ligand residue name"
    )
    parser.add_argument(
        "--benchmark_ligands", nargs="*", default=[], help="Benchmark ligand titles"
    )
    parser.add_argument("--outdir", required=True, help="Output reports directory")
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dock_dir = Path(args.docking_dir)
    analysis_dir = Path(args.analysis_dir) if args.analysis_dir else dock_dir

    if not args.plan:
        raise ValueError(
            "FATAL: --plan argument is required. Execution manifest must be passed explicitly."
        )

    plan_path = Path(args.plan).resolve()
    if not plan_path.exists():
        raise FileNotFoundError(
            f"FATAL: Execution plan JSON was not found at canonical path: {plan_path}"
        )

    receptors_ordered = []
    with open(plan_path, "r", encoding="utf-8") as pf:
        plan_data = json.load(pf)
        for rec in plan_data.get("receptors", []):
            receptors_ordered.append(rec["name"])
        for task in plan_data.get("tasks", []):
            t_name = task.get("name")
            if t_name and t_name not in receptors_ordered:
                receptors_ordered.append(t_name)

    if not receptors_ordered:
        raise ValueError(
            f"FATAL: Manifest {plan_path} does not specify any target receptors or tasks."
        )

    rows = []
    for rec_id in receptors_ordered:
        g_file = dock_dir / f"{rec_id}_dock.csv"
        m_file = analysis_dir / f"{rec_id}_mmgbsa-out.csv"
        g_resc_file = dock_dir / f"{rec_id}_glide_rescored.sdf"
        u_resc_file = (
            dock_dir / f"{rec_id}_unidock_rescored.sdf"
            if (dock_dir / f"{rec_id}_unidock_rescored.sdf").exists()
            else dock_dir / f"{rec_id}_unidock_dock.sdf"
        )

        glide_data = parse_glide_csv(g_file) if g_file.exists() else {}
        prime_data = parse_prime_csv(m_file) if m_file.exists() else {}
        glide_cnn_data = (
            parse_gnina_sdf(g_resc_file, prefix="Glide_") if g_resc_file.exists() else {}
        )
        unidock_data = (
            parse_gnina_sdf(u_resc_file, prefix="UniDock_") if u_resc_file.exists() else {}
        )

        titles = sorted(
            list(
                set(glide_data.keys())
                | set(prime_data.keys())
                | set(glide_cnn_data.keys())
                | set(unidock_data.keys())
            )
        )

        for t in titles:
            r = {"Title": t, "Receptor": rec_id, "Frame": rec_id}
            r.update(glide_data.get(t, {}))
            r.update(glide_cnn_data.get(t, {}))
            r.update(unidock_data.get(t, {}))
            r.update(prime_data.get(t, {}))
            rows.append(r)

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(
            "FATAL: No records found. Docking and/or rescoring failed entirely."
        )

    # Compute cross-engine Pose_RMSD first so it participates in the 3-axis consensus score
    df = compute_pose_consensus_rmsd(df, dock_dir)

    # Export full raw multi-frame metrics
    df_raw = df.copy()
    df_raw.to_csv(out_dir / "docking_combined_raw.csv", index=False)

    # Calculate 3-axis Orthogonal Consensus Score (Physics + CNNscore Validity + Pose_RMSD)
    df = calculate_consensus_score(df)

    # Pivot matrix: Ultimate_Score across all frames
    if "Ultimate_Score" in df.columns and "Frame" in df.columns:
        matrix_df = df.pivot_table(
            index="Title", columns="Frame", values="Ultimate_Score", aggfunc="max"
        )
        matrix_df.to_csv(out_dir / "cross_frame_matrix.csv")

    df_resolved = (
        df.sort_values("Ultimate_Score", ascending=False)
        .drop_duplicates(subset=["Title", "Frame"])
        .copy()
    )

    # Read separate Glide and Uni-Dock in-situ RMSD directly from calculated CSVs
    glide_rmsd_map = {}
    unidock_rmsd_map = {}

    # Deterministic lookup per target receptor from manifest without filesystem glob scanning
    for rec_id in receptors_ordered:
        for engine, target_map, col_name in [
            ("glide", glide_rmsd_map, "InSitu_RMSD_Glide_A"),
            ("unidock", unidock_rmsd_map, "InSitu_RMSD_UniDock_A"),
        ]:
            rmsd_csv = dock_dir / f"{rec_id}_{engine}_rmsd.csv"
            if rmsd_csv.exists() and rmsd_csv.stat().st_size > 0:
                r_df = pd.read_csv(rmsd_csv)
                for _, r_row in r_df.iterrows():
                    r_rec = str(r_row.get("Receptor", "")).strip()
                    r_title = str(r_row.get("Title", "")).strip()
                    if col_name in r_row and not pd.isna(r_row[col_name]):
                        target_map[(r_rec, r_title)] = float(r_row[col_name])

    glide_rmsds = []
    unidock_rmsds = []
    combined_rmsds = []

    for _, row in df_resolved.iterrows():
        rec_id = str(row["Receptor"]).strip()
        title = str(row["Title"]).strip()
        val_glide = glide_rmsd_map.get((rec_id, title), np.nan)
        val_unidock = unidock_rmsd_map.get((rec_id, title), np.nan)
        glide_rmsds.append(val_glide)
        unidock_rmsds.append(val_unidock)
        combined_rmsds.append(val_glide if not np.isnan(val_glide) else val_unidock)

    df_resolved["InSitu_RMSD_Glide_A"] = glide_rmsds
    df_resolved["InSitu_RMSD_UniDock_A"] = unidock_rmsds
    df_resolved["InSitu_RMSD_A"] = combined_rmsds

    # Delta score relative to native ligand if present
    native_targets = []
    if args.native_ligand:
        native_targets.append(args.native_ligand.lower())

    native_mask = df_resolved["Title"].str.lower().isin(native_targets)
    if native_targets and native_mask.any():
        native_score = df_resolved.loc[native_mask, "Ultimate_Score"].iloc[0]
        df_resolved["Delta_to_Native"] = df_resolved["Ultimate_Score"] - native_score

    df_resolved.reset_index(drop=True, inplace=True)
    df_resolved.index += 1
    df_resolved.index.name = "Rank"

    # Save master tables without numbering prefixes (both semantic and ranking aliases)
    df_resolved.to_csv(out_dir / "consensus_ranking.csv")
    df_resolved.to_csv(out_dir / "ranking.csv")

    # Save summary JSON
    with open(out_dir / "summary.json", "w") as fj:
        json.dump(df_resolved.reset_index().to_dict(orient="records"), fj, indent=2)

    print(
        f"Master consensus ranking generated with {len(df_resolved)} compounds in {out_dir / 'consensus_ranking.csv'}"
    )
    print(f"Cross-frame matrix saved in {out_dir / 'cross_frame_matrix.csv'}")
    print(
        f"Raw un-deduplicated metrics saved in {out_dir / 'docking_combined_raw.csv'}"
    )


if __name__ == "__main__":
    main()
