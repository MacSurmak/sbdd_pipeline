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

        # Search for glide and gnina poses in dock_dir with exact frame boundary
        patt = rf"_f{frame_id}(?:[^0-9]|$)"
        glide_sdfs = [
            f for f in dock_dir.glob("**/*dock*.sdf") if re.search(patt, f.name)
        ]
        gnina_sdfs = [
            f for f in dock_dir.glob("**/*gnina*.sdf") if re.search(patt, f.name)
        ]

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
        df["PR_CNN_Resc"] = (
            df["Resc_CNNaffinity"].rank(ascending=True, pct=True).fillna(0.0)
        )
    if "Dock_CNNaffinity" in df.columns:
        df["PR_GNINA_Dock"] = (
            df["Dock_CNNaffinity"].rank(ascending=True, pct=True).fillna(0.0)
        )

    pr_cols = [c for c in ["PR_Glide", "PR_MMGBSA", "PR_CNN_Resc"] if c in df.columns]
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
        "--docking_dir", help="Directory containing docking poses and tables"
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

    dock_dir = Path(args.docking_dir) if args.docking_dir else Path(".")
    analysis_dir = Path(args.analysis_dir) if args.analysis_dir else dock_dir

    # Find frames by searching CSVs and SDFs
    glide_csvs = list(dock_dir.glob("**/*glide*f*.csv")) or list(
        dock_dir.glob("**/*dock*.csv")
    )
    mmgbsa_csvs = list(analysis_dir.glob("**/*mmgbsa*.csv"))
    resc_sdfs = list(dock_dir.glob("**/*rescored*.sdf")) or list(
        dock_dir.glob("**/*gnina_score_only*.sdf")
    )
    dock_sdfs = list(dock_dir.glob("**/*gnina*dock*.sdf")) or list(
        dock_dir.glob("**/*docked*.sdf")
    )

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

        # Match files for frame with exact regex boundary (avoid _f1 matching _f10)
        pattern = rf"_f{fid}(?:[^0-9]|$)"
        g_files = [f for f in glide_csvs if re.search(pattern, f.name)]
        m_files = [f for f in mmgbsa_csvs if re.search(pattern, f.name)]
        r_files = [f for f in resc_sdfs if re.search(pattern, f.name)]
        d_files = [f for f in dock_sdfs if re.search(pattern, f.name)]

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

        titles = set(
            list(glide_data.keys())
            + list(prime_data.keys())
            + list(resc_data.keys())
            + list(dock_data.keys())
        )
        for t in titles:
            r = {"Title": t, "Frame": frame_tag}
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

    # Calculate in-situ redocking RMSD against reference pose if available
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign

        insitu_rmsds = []
        for _, row in df_resolved.iterrows():
            fid = str(row["Frame"]).replace("F", "").replace("f", "")
            title = row["Title"]
            ref_sdf = Path(f"prepared/{fid}_ref_native_prepared.sdf")
            if not ref_sdf.exists():
                ref_sdf = next(Path("prepared").glob(f"*{fid}*pose*.sdf"), None) or next(Path("prepared").glob(f"*{fid}*pose*.pdb"), None)

            rmsd_val = np.nan
            if ref_sdf and ref_sdf.exists():
                # Extract docked pose SDF block
                dock_candidates = list(dock_dir.glob(f"**/*f{fid}*dock*.sdf"))
                if dock_candidates:
                    with open(dock_candidates[0], "r", errors="ignore") as f:
                        for b in f.read().split("$$$$"):
                            lines = b.strip().splitlines()
                            if lines and (lines[0].strip() == title or title in lines[0]):
                                m_dock = Chem.RemoveHs(Chem.MolFromMolBlock(b.strip() + "\n$$$$\n"))
                                m_ref = Chem.RemoveHs(Chem.MolFromMolFile(str(ref_sdf)))
                                if m_dock and m_ref and m_dock.GetNumAtoms() == m_ref.GetNumAtoms():
                                    rmsd_val = round(rdMolAlign.GetBestRMS(m_dock, m_ref), 2)
                                break
            insitu_rmsds.append(rmsd_val)
        df_resolved["InSitu_RMSD_A"] = insitu_rmsds
    except Exception:
        df_resolved["InSitu_RMSD_A"] = np.nan

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
