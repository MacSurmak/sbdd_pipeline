#!/usr/bin/env python3
"""
bin/build_funnel_whitelist.py
Aggregates Glide and GNINA Rescore data to triage candidates using consensus (Physics + AI),
matching the exact logic from the original SBDD pipeline.
Outputs funnel_whitelist.csv with [Title, Frame].
"""
import argparse
import csv
import sys
from pathlib import Path
import pandas as pd


def parse_glide_csv(csv_path: Path) -> dict:
    data = {}
    if not csv_path.exists():
        return data
    with open(csv_path, "r", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("title") or row.get("Title") or row.get("s_m_title")
            if not title or title in data:
                continue

            score = None
            for k in ["r_i_docking_score", "docking_score"]:
                if k in row and row[k]:
                    try:
                        score = float(row[k])
                        break
                    except ValueError:
                        pass
            if score is None:
                for k in ["r_i_glide_gscore", "gscore"]:
                    if k in row and row[k]:
                        try:
                            score = float(row[k])
                            break
                        except ValueError:
                            pass
            if score is not None:
                data[title] = {"GlideScore": score}
    return data


def parse_gnina_sdf(sdf_path: Path, prefix: str = "") -> dict:
    import re
    data = {}
    if not sdf_path or not sdf_path.exists():
        return data
    with open(sdf_path, "r", errors="ignore") as f:
        content = f.read()
    for entry in content.split("$$$$"):
        lines = entry.strip().splitlines()
        if not lines:
            continue

        raw_title = lines[0].strip()
        prop_title = None
        props = {}

        for i, line in enumerate(lines):
            line_str = line.strip()
            if line_str.startswith("> <") and line_str.endswith(">"):
                tag_clean = line_str[3:-1].strip().lower()
                if i + 1 < len(lines):
                    val_str = lines[i + 1].strip()
                    if tag_clean in ["_name", "s_m_title", "title", "compound_id"]:
                        prop_title = val_str
                    else:
                        try:
                            val = float(val_str)
                            if tag_clean == "cnnscore":
                                props[f"{prefix}CNNscore"] = val
                            elif tag_clean == "cnnaffinity":
                                props[f"{prefix}CNNaffinity"] = val
                            elif tag_clean in ["minimizedaffinity", "affinity"]:
                                props[f"{prefix}VinaAff"] = val
                        except ValueError:
                            pass

        title = prop_title or raw_title
        if not title:
            continue

        title_clean = re.sub(r"(_model_\d+|_dock_\d+|_\d+)$", "", title).strip()

        if props:
            if title not in data:
                data[title] = props
            if title_clean and title_clean not in data:
                data[title_clean] = props
    return data


def main():
    parser = argparse.ArgumentParser(description="Build screening funnel whitelist")
    parser.add_argument("--glide_csvs", nargs="+", default=[], help="Glide docking CSV files")
    parser.add_argument("--gnina_sdfs", nargs="+", default=[], help="GNINA rescore SDF files")
    parser.add_argument("--top_k_frames", type=int, default=3, help="Max frames per ligand")
    parser.add_argument("--percentile_cutoff", type=float, default=0.50, help="Fraction of top ligands (e.g. 0.50)")
    parser.add_argument("--output", required=True, help="Output whitelist CSV path")
    args = parser.parse_args()

    # Map frame id to files
    import re
    glide_by_frame = {}
    for f in args.glide_csvs:
        p = Path(f)
        # expected format: ..._f<frame>... or ..._f<frame>.csv
        stem = p.stem
        parts = stem.split("_f")
        frame_id = parts[-1].split("_")[0] if len(parts) > 1 else "1"
        m = re.search(r"_f(\d+)", p.stem)
        frame_id = m.group(1) if m else "1"
        glide_by_frame[frame_id] = p

    gnina_by_frame = {}
    for f in args.gnina_sdfs:
        p = Path(f)
        parts = p.stem.split("_f")
        frame_id = parts[-1].split("_")[0] if len(parts) > 1 else "1"
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
                "Resc_CNNaffinity": resc_data.get(title, {}).get("Resc_CNNaffinity", None)
            }
            all_data.append(entry)

    df = pd.DataFrame(all_data)
    out_csv = Path(args.output)
    if df.empty:
        print("Warning: No docking data available to build funnel whitelist. Outputting empty whitelist.")
        pd.DataFrame(columns=["Title", "Frame"]).to_csv(out_csv, index=False)
        sys.exit(0)

    # 1. Percentile Ranks (Glide lower is better, CNN higher is better)
    df["PR_Glide"] = df["GlideScore"].rank(ascending=False, pct=True)
    if df["Resc_CNNaffinity"].notna().any():
        df["PR_CNN"] = df["Resc_CNNaffinity"].rank(ascending=True, pct=True).fillna(df["PR_Glide"])
        df["Triage_Score"] = (df["PR_Glide"] + df["PR_CNN"]) / 2.0
    else:
        df["Triage_Score"] = df["PR_Glide"]

    # 2. Cutoff top percentile (e.g. top 50%)
    best_triage = df.groupby("Title")["Triage_Score"].max()
    threshold = best_triage.quantile(1.0 - args.percentile_cutoff)
    kept_ligands = set(best_triage[best_triage >= threshold].index)

    df_filtered = df[df["Title"].isin(kept_ligands)]

    # 3. Top K frames per ligand
    whitelist = df_filtered.sort_values("Triage_Score", ascending=False).groupby("Title").head(args.top_k_frames)

    whitelist[["Title", "Frame"]].to_csv(out_csv, index=False)
    print(f"Screening Funnel: Retained {len(kept_ligands)}/{len(best_triage)} ligands across {len(whitelist)} tasks -> {out_csv}")


if __name__ == "__main__":
    main()
