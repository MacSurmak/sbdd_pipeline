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
import csv
import re
from pathlib import Path
import pandas as pd
import numpy as np


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


def parse_prime_csv(csv_path: Path) -> dict:
    data = {}
    if not csv_path.exists():
        return data
    with open(csv_path, "r", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("title") or row.get("Title") or row.get("s_m_title")
            if not title or title in data:
                continue
            props = {}
            for k, v in row.items():
                if not k or not v:
                    continue
                k_clean = k.strip()
                if k_clean in ["r_psp_MMGBSA_dG_Bind", "MMGBSA dG Bind"]:
                    try:
                        props["MMGBSA"] = float(v)
                        break
                    except ValueError:
                        pass

            for k, v in row.items():
                if not k or not v:
                    continue
                if k.strip() in ["r_psp_Lig_Strain_Energy", "Lig Strain Energy"]:
                    try:
                        props["Lig_Strain"] = float(v)
                        break
                    except ValueError:
                        pass

            if props:
                data[title] = props
    return data


def parse_gnina_sdf(sdf_path: Path, prefix: str = "") -> dict:
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


def compute_pose_consensus_rmsd(df: pd.DataFrame, dock_dir: Path):
    """Calculates in-place heavy-atom RMSD between Glide and independent GNINA dock poses."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
    except ImportError:
        df["Pose_RMSD"] = np.nan
        return df

    def _extract_sdf_block(filepath: Path, target_title: str) -> str:
        if not filepath.exists():
            return ""
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
        for block in content.split("$$$$"):
            lines = block.strip().splitlines()
            if not lines:
                continue
            title = lines[0].strip()
            if title == target_title:
                return block.strip() + "\n$$$$\n"
        return ""

    rmsd_list = []
    for _, row in df.iterrows():
        frame_id = str(row["Frame"]).replace("F", "").replace("f", "")
        title = row["Title"]

        # Search for glide and gnina poses in dock_dir
        glide_sdfs = list(dock_dir.glob(f"**/*f{frame_id}*dock*.sdf"))
        gnina_sdfs = list(dock_dir.glob(f"**/*f{frame_id}*gnina*.sdf"))
        # Search for glide and gnina poses in dock_dir with exact frame boundary
        patt = rf"_f{frame_id}(?:[^0-9]|$)"
        glide_sdfs = [f for f in dock_dir.glob("**/*dock*.sdf") if re.search(patt, f.name)]
        gnina_sdfs = [f for f in dock_dir.glob("**/*gnina*.sdf") if re.search(patt, f.name)]

        glide_sdf = glide_sdfs[0] if glide_sdfs else None
        gnina_sdf = gnina_sdfs[0] if gnina_sdfs else None

        if glide_sdf and gnina_sdf:
            block_g = _extract_sdf_block(glide_sdf, title)
            block_n = _extract_sdf_block(gnina_sdf, title)
            if block_g and block_n:
                try:
                    mg = Chem.RemoveHs(Chem.MolFromMolBlock(block_g))
                    mn = Chem.RemoveHs(Chem.MolFromMolBlock(block_n))
                    if mg and mn and mg.GetNumAtoms() == mn.GetNumAtoms():
                        val = rdMolAlign.GetBestRMS(mg, mn)
                        rmsd_list.append(round(val, 2))
                        continue
                except Exception:
                    pass
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
        df["PR_CNN_Resc"] = df["Resc_CNNaffinity"].rank(ascending=True, pct=True).fillna(0.0)
    if "Dock_CNNaffinity" in df.columns:
        df["PR_GNINA_Dock"] = df["Dock_CNNaffinity"].rank(ascending=True, pct=True).fillna(0.0)

    pr_cols = [c for c in ["PR_Glide", "PR_MMGBSA", "PR_CNN_Resc"] if c in df.columns]
    if not pr_cols:
        pr_cols = [c for c in ["PR_GNINA_Dock"] if c in df.columns]

    if not pr_cols:
        df["Ultimate_Score"] = 0.0
    else:
        df["Ultimate_Score"] = (df[pr_cols].sum(axis=1) / len(pr_cols)) * 100.0
    return df


def main():
    parser = argparse.ArgumentParser(description="Aggregate project results into master ranking")
    parser.add_argument("--docking_dir", help="Directory containing docking poses and tables")
    parser.add_argument("--analysis_dir", help="Directory containing MM-GBSA csv files")
    parser.add_argument("--outdir", required=True, help="Output reports directory")
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dock_dir = Path(args.docking_dir) if args.docking_dir else Path(".")
    analysis_dir = Path(args.analysis_dir) if args.analysis_dir else dock_dir

    # Find frames by searching CSVs and SDFs
    glide_csvs = list(dock_dir.glob("**/*glide*f*.csv")) or list(dock_dir.glob("**/*dock*.csv"))
    mmgbsa_csvs = list(analysis_dir.glob("**/*mmgbsa*.csv"))
    resc_sdfs = list(dock_dir.glob("**/*rescored*.sdf")) or list(dock_dir.glob("**/*gnina_score_only*.sdf"))
    dock_sdfs = list(dock_dir.glob("**/*gnina_dock*.sdf"))

    # Map by frame ID
    all_frames = set()
    for p in glide_csvs + mmgbsa_csvs + resc_sdfs + dock_sdfs:
        match = re.search(r"_f(\d+)", p.stem)
        if match:
            all_frames.add(match.group(1))

    rows = []
    if not all_frames:
        # Fallback to single frame 1
        all_frames = {"1"}

    for fid in sorted(all_frames, key=lambda x: int(x) if x.isdigit() else x):
        frame_tag = f"F{fid}"

        # Match files for frame
        g_file = next((f for f in glide_csvs if f"_f{fid}" in f.name), None)
        m_file = next((f for f in mmgbsa_csvs if f"_f{fid}" in f.name), None)
        r_file = next((f for f in resc_sdfs if f"_f{fid}" in f.name), None)
        d_file = next((f for f in dock_sdfs if f"_f{fid}" in f.name), None)
        # Match files for frame with exact regex boundary (avoid _f1 matching _f10)
        pattern = rf"_f{fid}(?:[^0-9]|$)"
        g_files = [f for f in glide_csvs if re.search(pattern, f.name)]
        m_files = [f for f in mmgbsa_csvs if re.search(pattern, f.name)]
        r_files = [f for f in resc_sdfs if re.search(pattern, f.name)]
        d_files = [f for f in dock_sdfs if re.search(pattern, f.name)]

        glide_data = parse_glide_csv(g_file) if g_file else {}
        prime_data = parse_prime_csv(m_file) if m_file else {}
        resc_data = parse_gnina_sdf(r_file, prefix="Resc_") if r_file else {}
        dock_data = parse_gnina_sdf(d_file, prefix="Dock_") if d_file else {}
        glide_data = {}
        for f in g_files:
            glide_data.update(parse_glide_csv(f))

        prime_data = {}
        for f in m_files:
            prime_data.update(parse_prime_csv(f))

        resc_data = {}
        for f in r_files:
            resc_data.update(parse_gnina_sdf(f, prefix="Resc_"))

        dock_data = {}
        for f in d_files:
            dock_data.update(parse_gnina_sdf(f, prefix="Dock_"))

        titles = set(list(glide_data.keys()) + list(prime_data.keys()) + list(resc_data.keys()) + list(dock_data.keys()))
        for t in titles:
            r = {"Title": t, "Frame": frame_tag}
            r.update(glide_data.get(t, {}))
            r.update(prime_data.get(t, {}))
            r.update(resc_data.get(t, {}))
            r.update(dock_data.get(t, {}))
            rows.append(r)

    df = pd.DataFrame(rows)
    if df.empty:
        print("Warning: No records found. Generating empty templates.")
        pd.DataFrame(columns=["Rank", "Title", "Frame", "Ultimate_Score"]).to_csv(out_dir / "ranking.csv", index=False)
        pd.DataFrame(columns=["Title"]).to_csv(out_dir / "cross_frame_matrix.csv", index=False)
        with open(out_dir / "summary.json", "w") as fj:
            json.dump([], fj)
        sys.exit(0)

    # Calculate Consensus Score
    df = calculate_consensus_score(df)

    # Pivot matrix: Ultimate_Score across all frames
    if "Ultimate_Score" in df.columns and "Frame" in df.columns:
        matrix_df = df.pivot_table(index="Title", columns="Frame", values="Ultimate_Score", aggfunc="max")
        matrix_df.to_csv(out_dir / "cross_frame_matrix.csv")

    # Resolved top hit per ligand
    df_resolved = df.sort_values("Ultimate_Score", ascending=False).drop_duplicates(subset="Title").copy()
    df_resolved = compute_pose_consensus_rmsd(df_resolved, dock_dir)

    df_resolved.reset_index(drop=True, inplace=True)
    df_resolved.index += 1
    df_resolved.index.name = "Rank"

    # Save ranking.csv
    df_resolved.to_csv(out_dir / "ranking.csv")

    # Save summary.json
    with open(out_dir / "summary.json", "w") as fj:
        json.dump(df_resolved.reset_index().to_dict(orient="records"), fj, indent=2)

    print(f"Master consensus ranking generated with {len(df_resolved)} compounds in {out_dir / 'ranking.csv'}")
    print(f"Cross-frame matrix saved in {out_dir / 'cross_frame_matrix.csv'}")


if __name__ == "__main__":
    main()
