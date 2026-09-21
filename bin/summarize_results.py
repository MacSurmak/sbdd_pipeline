#!/usr/bin/env python3
"""
File: bin/summarize_results.py
Description: Aggregates raw docking, MM-GBSA, and MD metrics into a unified ranking.csv and summary.json.
"""
import argparse
import sys
import json
import csv
from pathlib import Path
import pandas as pd

def parse_gnina_sdf(sdf_path: Path):
    """Extracts scores from GNINA SDF poses."""
    poses = []
    current_title = ""
    current_props = {}
    
    with open(sdf_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            if i == 0 or (i > 0 and lines[i-1].strip() == "$$$$"):
                current_title = line.strip()
            if line.startswith("> <"):
                key = line.strip()[3:-1]
                val = lines[i+1].strip() if i+1 < len(lines) else ""
                current_props[key] = val
                i += 1
            elif line.strip() == "$$$$":
                if current_title:
                    poses.append({
                        "Title": current_title,
                        "CNNscore": float(current_props.get("CNNscore", 0.0)) if current_props.get("CNNscore") else None,
                        "CNNaffinity": float(current_props.get("CNNaffinity", 0.0)) if current_props.get("CNNaffinity") else None,
                        "minimizedAffinity": float(current_props.get("minimizedAffinity", 0.0)) if current_props.get("minimizedAffinity") else None
                    })
                cand_title = (
                    current_title 
                    or current_props.get("Title") 
                    or current_props.get("_Name") 
                    or current_props.get("Name") 
                    or f"Compound_{len(poses) + 1}"
                )
                poses.append({
                    "Title": cand_title,
                    "CNNscore": float(current_props.get("CNNscore", 0.0)) if current_props.get("CNNscore") else None,
                    "CNNaffinity": float(current_props.get("CNNaffinity", 0.0)) if current_props.get("CNNaffinity") else None,
                    "minimizedAffinity": float(current_props.get("minimizedAffinity", 0.0)) if current_props.get("minimizedAffinity") else None
                })
                current_title = ""
                current_props = {}
            i += 1
            
    df = pd.DataFrame(poses)
    if not df.empty:
        # Sort by CNNaffinity descending (higher affinity = better) or minimizedAffinity ascending
        sort_col = "CNNaffinity" if "CNNaffinity" in df.columns and df["CNNaffinity"].notna().any() else "minimizedAffinity"
        ascending = False if sort_col == "CNNaffinity" else True
        df = df.sort_values(sort_col, ascending=ascending)
    return df

def main():
    parser = argparse.ArgumentParser(description="Aggregate project results into master ranking")
    parser.add_argument("--docking_dir", help="Directory containing docking poses and tables")
    parser.add_argument("--analysis_dir", help="Directory containing MM-GBSA csv files")
    parser.add_argument("--md_dir", help="Directory containing MD results")
    parser.add_argument("--outdir", required=True, help="Output reports directory")
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    master_rows = {}

    # 1. Ingest Docking results
    if args.docking_dir and Path(args.docking_dir).exists():
        dock_p = Path(args.docking_dir)
        # Parse CSV tables if Glide
        for csv_f in dock_p.glob("**/*.csv"):
            try:
                df = pd.read_csv(csv_f)
                for _, r in df.iterrows():
                    title = str(r.get("Title") or r.get("s_m_title") or "").strip()
                    if not title: continue
                    if title not in master_rows: master_rows[title] = {"Title": title}
                    if "r_i_glide_gscore" in r: master_rows[title]["GlideScore"] = float(r["r_i_glide_gscore"])
                    if "r_i_glide_docking_score" in r: master_rows[title]["DockingScore"] = float(r["r_i_glide_docking_score"])
                    if "RMSD_Heavy" in r: master_rows[title]["RMSD_A"] = float(r["RMSD_Heavy"])
                    for k in ["r_psp_MMGBSA_dG_Bind", "MMGBSA_dG_Bind", "dG_Bind"]:
                        if k in r and pd.notna(r[k]):
                            master_rows[title]["Prime_dG_Bind_kcal"] = float(r[k])
                            break
            except Exception:
                pass

        # Parse SDFs if GNINA
        for sdf_f in dock_p.glob("**/*.sdf"):
            try:
                df_gnina = parse_gnina_sdf(sdf_f)
                for _, r in df_gnina.iterrows():
                    title = r["Title"]
                    if not title: continue
                    if title not in master_rows: master_rows[title] = {"Title": title}
                    if r["CNNaffinity"] is not None: master_rows[title]["CNNaffinity"] = r["CNNaffinity"]
                    if r["CNNscore"] is not None: master_rows[title]["CNNscore"] = r["CNNscore"]
                    if r["minimizedAffinity"] is not None: master_rows[title]["GNINA_Affinity"] = r["minimizedAffinity"]
            except Exception:
                pass

    # 2. Ingest Prime MM-GBSA
    if args.analysis_dir and Path(args.analysis_dir).exists():
        for csv_f in Path(args.analysis_dir).glob("**/*mmgbsa*.csv"):
            try:
                df = pd.read_csv(csv_f)
                for _, r in df.iterrows():
                    title = str(r.get("Title") or r.get("s_m_title") or "").strip()
                    if not title: continue
                    if title not in master_rows: master_rows[title] = {"Title": title}
                    for k in ["r_psp_MMGBSA_dG_Bind", "MMGBSA_dG_Bind", "dG_Bind"]:
                        if k in r and pd.notna(r[k]):
                            master_rows[title]["Prime_dG_Bind_kcal"] = float(r[k])
                            break
            except Exception:
                pass

    # 3. Ingest MD Metrics
    if args.md_dir and Path(args.md_dir).exists():
        for json_f in Path(args.md_dir).glob("**/md_metrics_summary.json"):
            try:
                target_name = json_f.parent.name
                with open(json_f) as jf:
                    data = json.load(jf)
                matched_key = next((k for k in master_rows if k.lower() == target_name.lower()), target_name)
                if matched_key not in master_rows: master_rows[matched_key] = {"Title": matched_key}
                master_rows[matched_key]["MD_Prot_RMSD_Mean"] = data.get("Prot_BB_RMSD_Mean_A")
                master_rows[matched_key]["MD_Lig_RMSD_Mean"] = data.get("Lig_Heavy_RMSD_Mean_A")
            except Exception:
                pass

    df_out = pd.DataFrame(list(master_rows.values()))
    if not df_out.empty:
        # Determine ranking sort column
        if "Prime_dG_Bind_kcal" in df_out.columns:
            df_out = df_out.sort_values("Prime_dG_Bind_kcal", ascending=True)
        elif "CNNaffinity" in df_out.columns:
            df_out = df_out.sort_values("CNNaffinity", ascending=False)
        elif "GlideScore" in df_out.columns:
            df_out = df_out.sort_values("GlideScore", ascending=True)

        df_out.reset_index(drop=True, inplace=True)
        df_out.index += 1
        df_out.index.name = "Rank"
        df_out.to_csv(out_dir / "ranking.csv")
        
        with open(out_dir / "summary.json", "w") as fj:
            json.dump(df_out.reset_index().to_dict(orient="records"), fj, indent=2)
        print(f"Master ranking generated with {len(df_out)} candidates in {out_dir / 'ranking.csv'}")
    else:
        print("Warning: No records found to summarize.")
        pd.DataFrame(columns=["Rank", "Title", "Score"]).to_csv(out_dir / "ranking.csv", index=False)
        with open(out_dir / "summary.json", "w") as fj:
            json.dump([], fj)
        print("Warning: No records found to summarize. Empty template generated.")

if __name__ == "__main__":
    main()

