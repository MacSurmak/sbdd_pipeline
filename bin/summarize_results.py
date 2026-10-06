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
import sys
import json
import re
from pathlib import Path
import pandas as pd
import numpy as np

# Add script directory to sys.path for absolute imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.parsers import parse_glide_csv, parse_prime_csv, parse_gnina_sdf


def compute_pose_consensus_rmsd(df: pd.DataFrame, dock_dir: Path):
    """Calculates in-place heavy-atom RMSD between Glide and GNINA dock poses using explicit file paths."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
    except ImportError:
        df["Pose_RMSD"] = np.nan
        return df

    def _extract_mol(filepath: Path, target_title: str):
        if not filepath.exists() or filepath.stat().st_size == 0:
            return None
        suppl = Chem.SDMolSupplier(str(filepath), removeHs=True, sanitize=False)
        for m in suppl:
            if m is not None:
                m_title = m.GetProp("_Name").strip() if m.HasProp("_Name") else ""
                if m_title == target_title:
                    return m
        return None

    rmsd_list = []
    for _, row in df.iterrows():
        rec_id = str(row["Receptor"]).strip()
        title = str(row["Title"]).strip()

        glide_sdf = dock_dir / f"{rec_id}_dock.sdf"
        gnina_sdf = dock_dir / f"{rec_id}_gnina_docked_all.sdf"

        if glide_sdf.exists() and gnina_sdf.exists():
            mg = _extract_mol(glide_sdf, title)
            mn = _extract_mol(gnina_sdf, title)
            if (
                mg is not None
                and mn is not None
                and mg.GetNumAtoms() == mn.GetNumAtoms()
            ):
                val = rdMolAlign.GetBestRMS(mg, mn)
                rmsd_list.append(round(val, 2))
                continue
        rmsd_list.append(np.nan)

    df["Pose_RMSD"] = rmsd_list
    return df


def calculate_consensus_score(df: pd.DataFrame) -> pd.DataFrame:
    """Computes Percentile Rank consensus score (0-100) exactly matching original pipeline."""
    if "GlideScore" in df.columns:
        df["PR_Glide"] = df["GlideScore"].rank(ascending=False, pct=True).fillna(0.0)
    if "MMGBSA" in df.columns:
        df["PR_MMGBSA"] = df["MMGBSA"].rank(ascending=False, pct=True).fillna(0.0)
    if "Resc_CNNaffinity" in df.columns:
        df["PR_CNN_Resc"] = (
            df["Resc_CNNaffinity"].rank(ascending=True, pct=True).fillna(0.0)
        )
    if "Dock_CNNaffinity" in df.columns:
        df["PR_GNINA_Dock"] = (
            df["Dock_CNNaffinity"].rank(ascending=True, pct=True).fillna(0.0)
        )

    # Gracefully adapt consensus ranking whether MM-GBSA was executed or deferred
    candidate_cols = ["PR_Glide", "PR_MMGBSA", "PR_CNN_Resc"]
    pr_cols = [c for c in candidate_cols if c in df.columns and df[c].notna().any()]
    if not pr_cols:
        pr_cols = [c for c in ["PR_GNINA_Dock"] if c in df.columns]

    if not pr_cols:
        df["Ultimate_Score"] = 0.0
    else:
        df["Ultimate_Score"] = (df[pr_cols].sum(axis=1) / len(pr_cols)) * 100.0
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
        r_file = dock_dir / f"{rec_id}_rescored.sdf"
        d_file = dock_dir / f"{rec_id}_gnina_docked_all.sdf"

        glide_data = parse_glide_csv(g_file) if g_file.exists() else {}
        prime_data = parse_prime_csv(m_file) if m_file.exists() else {}
        resc_data = parse_gnina_sdf(r_file, prefix="Resc_") if r_file.exists() else {}
        dock_data = parse_gnina_sdf(d_file, prefix="Dock_") if d_file.exists() else {}

        titles = sorted(
            list(
                set(glide_data.keys())
                | set(prime_data.keys())
                | set(resc_data.keys())
                | set(dock_data.keys())
            )
        )

        for t in titles:
            r = {"Title": t, "Receptor": rec_id, "Frame": rec_id}
            r.update(glide_data.get(t, {}))
            r.update(prime_data.get(t, {}))
            r.update(resc_data.get(t, {}))
            r.update(dock_data.get(t, {}))
            rows.append(r)

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(
            "FATAL: No records found. Docking and/or rescoring failed entirely."
        )

    # Export full raw multi-frame metrics before deduplication
    df_raw = df.copy()
    df_raw.to_csv(out_dir / "docking_combined_raw.csv", index=False)

    # Calculate Consensus Score
    df = calculate_consensus_score(df)

    # Pivot matrix: Ultimate_Score across all frames
    if "Ultimate_Score" in df.columns and "Frame" in df.columns:
        matrix_df = df.pivot_table(
            index="Title", columns="Frame", values="Ultimate_Score", aggfunc="max"
        )
        matrix_df.to_csv(out_dir / "cross_frame_matrix.csv")

    # Resolved top hit per ligand (best conformational frame)
    df_resolved = (
        df.sort_values("Ultimate_Score", ascending=False)
        .drop_duplicates(subset=["Title", "Frame"])
        .copy()
    )
    df_resolved = compute_pose_consensus_rmsd(df_resolved, dock_dir)

    # Read separate Glide and GNINA in-situ RMSD directly from calculated CSVs
    glide_rmsd_map = {}
    gnina_rmsd_map = {}

    # Deterministic lookup per target receptor from manifest without filesystem glob scanning
    for rec_id in receptors_ordered:
        for engine, target_map, col_name in [
            ("glide", glide_rmsd_map, "InSitu_RMSD_Glide_A"),
            ("gnina", gnina_rmsd_map, "InSitu_RMSD_GNINA_A"),
        ]:
            rmsd_csv = dock_dir / f"{rec_id}_{engine}_rmsd.csv"
            if rmsd_csv.exists() and rmsd_csv.stat().st_size > 0:
                try:
                    r_df = pd.read_csv(rmsd_csv)
                    for _, r_row in r_df.iterrows():
                        r_rec = str(r_row.get("Receptor", "")).strip()
                        r_title = str(r_row.get("Title", "")).strip()
                        if col_name in r_row and not pd.isna(r_row[col_name]):
                            target_map[(r_rec, r_title)] = float(r_row[col_name])
                except Exception as e:
                    print(
                        f"Warning: Failed to parse RMSD CSV {rmsd_csv}: {e}",
                        file=sys.stderr,
                    )

    glide_rmsds = []
    gnina_rmsds = []
    combined_rmsds = []

    for _, row in df_resolved.iterrows():
        rec_id = str(row["Receptor"]).strip()
        title = str(row["Title"]).strip()
        val_glide = glide_rmsd_map.get((rec_id, title), np.nan)
        val_gnina = gnina_rmsd_map.get((rec_id, title), np.nan)
        glide_rmsds.append(val_glide)
        gnina_rmsds.append(val_gnina)
        # Primary RMSD defaults to Glide if available, otherwise GNINA
        combined_rmsds.append(val_glide if not np.isnan(val_glide) else val_gnina)

    df_resolved["InSitu_RMSD_Glide_A"] = glide_rmsds
    df_resolved["InSitu_RMSD_GNINA_A"] = gnina_rmsds
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
