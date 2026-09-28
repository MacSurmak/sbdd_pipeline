#!/usr/bin/env python3
"""
bin/build_funnel_whitelist.py
Aggregates Glide and GNINA Rescore data to triage candidates using consensus (Physics + AI),
matching the exact logic from the original SBDD pipeline.
Outputs funnel_whitelist.csv with [Title, Frame].
"""
import argparse
import sys
from pathlib import Path
import pandas as pd

# Add script directory to sys.path for absolute imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.parsers import parse_glide_csv, parse_gnina_sdf


def main():
    parser = argparse.ArgumentParser(description="Build screening funnel whitelist")
    parser.add_argument(
        "--glide_csvs", nargs="+", default=[], help="Glide docking CSV files"
    )
    parser.add_argument(
        "--gnina_sdfs", nargs="+", default=[], help="GNINA rescore SDF files"
    )
    parser.add_argument(
        "--top_k_frames", type=int, default=3, help="Max frames per ligand"
    )
    parser.add_argument(
        "--percentile_cutoff",
        type=float,
        default=0.50,
        help="Fraction of top ligands (e.g. 0.50)",
    )
    parser.add_argument(
        "--benchmark_ligands",
        nargs="*",
        default=[],
        help="Benchmark/control ligands to always preserve",
    )
    parser.add_argument(
        "--native_ligand",
        default="",
        help="Native ligand residue name or title to always preserve",
    )
    parser.add_argument("--output", required=True, help="Output whitelist CSV path")
    args = parser.parse_args()

    # Map frame id to files
    import re

    glide_by_frame = {}
    for f in args.glide_csvs:
        p = Path(f)
        m = re.search(r"_f(\d+)", p.stem)
        frame_id = m.group(1) if m else "1"
        glide_by_frame[frame_id] = p

    gnina_by_frame = {}
    for f in args.gnina_sdfs:
        p = Path(f)
        m = re.search(r"_f(\d+)", p.stem)
        frame_id = m.group(1) if m else "1"
        gnina_by_frame[frame_id] = p

    all_data = []
    for frame_id, g_csv in glide_by_frame.items():
        glide_data = parse_glide_csv(g_csv)
        resc_sdf = gnina_by_frame.get(frame_id)
        resc_data = parse_gnina_sdf(resc_sdf, prefix="Resc_") if resc_sdf else {}

        for title, scores in glide_data.items():
            entry = {
                "Title": title,
                "Frame": frame_id,
                "GlideScore": scores["GlideScore"],
                "Resc_CNNaffinity": resc_data.get(title, {}).get(
                    "Resc_CNNaffinity", None
                ),
            }
            all_data.append(entry)

    df = pd.DataFrame(all_data)
    out_csv = Path(args.output)
    if df.empty:
        raise ValueError("FATAL: No docking data available to build funnel whitelist.")

    # 1. Percentile Ranks (Glide lower is better, CNN higher is better)
    df["PR_Glide"] = df["GlideScore"].rank(ascending=False, pct=True)
    if df["Resc_CNNaffinity"].notna().any():
        df["PR_CNN"] = (
            df["Resc_CNNaffinity"].rank(ascending=True, pct=True).fillna(df["PR_Glide"])
        )
        df["Triage_Score"] = (df["PR_Glide"] + df["PR_CNN"]) / 2.0
    else:
        df["Triage_Score"] = df["PR_Glide"]

    # 2. Cutoff top percentile (e.g. top 50%) and unconditionally keep benchmark controls and native ligand
    best_triage = df.groupby("Title")["Triage_Score"].max()
    threshold = best_triage.quantile(1.0 - args.percentile_cutoff)
    passing_score = set(best_triage[best_triage >= threshold].index)
    forced_set = set(args.benchmark_ligands or [])
    if args.native_ligand:
        forced_set.add(args.native_ligand)
        for t in df["Title"]:
            if args.native_ligand.lower() == t.lower():
                forced_set.add(t)

    kept_ligands = passing_score | (forced_set & set(df["Title"]))

    df_filtered = df[df["Title"].isin(kept_ligands)]

    # 3. Top K frames per ligand
    whitelist = (
        df_filtered.sort_values("Triage_Score", ascending=False)
        .groupby("Title")
        .head(args.top_k_frames)
    )

    whitelist[["Title", "Frame"]].to_csv(out_csv, index=False)
    print(
        f"Screening Funnel: Retained {len(kept_ligands)}/{len(best_triage)} ligands across {len(whitelist)} tasks -> {out_csv}"
    )


if __name__ == "__main__":
    main()
