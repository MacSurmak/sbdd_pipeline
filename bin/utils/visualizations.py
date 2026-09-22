"""
File: src/utils/visualizations.py
Description: Generation of consensus scatter plots and MPO heatmaps.
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

def generate_plots(df: pd.DataFrame, out_dir: Path, top_n: int):
    if df.empty: return
    sns.set_theme(style="whitegrid")
    
    # 1. Consensus Scatter Plot
    if all(c in df.columns for c in ["GlideScore", "MMGBSA", "Ultimate_Score"]):
        plt.figure(figsize=(10, 8))
        cnn_col = "Resc_CNNaffinity" if "Resc_CNNaffinity" in df.columns else "Dock_CNNaffinity" if "Dock_CNNaffinity" in df.columns else None
        sizes = df[cnn_col] * 10 if cnn_col else 50
        scatter = plt.scatter(
            x=df["GlideScore"], y=df["MMGBSA"], s=sizes, 
            c=df["Ultimate_Score"], cmap="viridis", alpha=0.8, edgecolor="k"
        )
        native_df = df[df["Is_Native"] == True]
        if not native_df.empty:
            plt.axvline(x=native_df["GlideScore"].iloc[0], color='r', linestyle='--', alpha=0.5)
            plt.axhline(y=native_df["MMGBSA"].iloc[0], color='r', linestyle='--', alpha=0.5)
            plt.scatter(
                native_df["GlideScore"], native_df["MMGBSA"], 
                s=150, facecolors='none', edgecolors='r', linewidth=2, label="Native Baseline"
            )
            plt.legend()
            
        plt.colorbar(scatter, label="Ultimate Score")
        plt.xlabel("GlideScore (kcal/mol)")
        plt.ylabel("Prime MM-GBSA (kcal/mol)")
        plt.title("Thermodynamic Landscape & Consensus Score")
        plt.tight_layout()
        plt.savefig(out_dir / "04_consensus_scatter.png", dpi=300)
        plt.close()

    # 2. MPO Heatmap
    pr_cols = [c for c in ["PR_Glide", "PR_MMGBSA", "PR_CNN_Resc", "PR_GNINA_Dock"] if c in df.columns]
    if pr_cols:
        top_df = df.head(top_n).set_index("Title")[pr_cols]
        plt.figure(figsize=(8, max(4, int(len(top_df) * 0.3))))
        sns.heatmap(top_df, cmap="RdYlGn", vmin=0, vmax=1, annot=True, fmt=".2f", cbar_kws={'label': 'Percentile Rank'})
        plt.title(f"Multiparameter Optimization (Top {len(top_df)})")
        plt.tight_layout()
        plt.savefig(out_dir / "03_mpo_heatmap.png", dpi=300)
        plt.close()

def plot_md_rmsd(df: pd.DataFrame, out_file: Path):
    if df.empty: return
    plt.figure(figsize=(10, 5))
    plt.plot(df["Time_ns"], df["Prot_BB_RMSD"], label="Protein Backbone", color="gray", alpha=0.8)
    plt.plot(df["Time_ns"], df["Lig_Heavy_RMSD"], label="Ligand Heavy", color="red", linewidth=2)
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("Molecular Dynamics RMSD")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

def plot_md_rmsf(df: pd.DataFrame, df_contacts: pd.DataFrame, out_file: Path):
    if df.empty: return
    plt.figure(figsize=(12, 4))
    plt.plot(df["ResID"], df["RMSF"], color="blue", linewidth=1.2)
    
    # Highlight high-frequency contact residues
    if not df_contacts.empty:
        occ_col = "Total_%" if "Total_%" in df_contacts.columns else "Occupancy_pct"
        if occ_col in df_contacts.columns:
            high_contacts = df_contacts[df_contacts[occ_col] >= 20.0]["ResID"]
            for resid in high_contacts:
                val = df[df["ResID"] == resid]["RMSF"].values
                if len(val) > 0:
                    plt.axvline(x=resid, color='orange', alpha=0.4, linestyle='--')
                
    plt.xlabel("Residue Index")
    plt.ylabel("RMSF (Å)")
    plt.title("Protein C-alpha RMSF (Interacting residues highlighted)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()

def plot_mmpbsa_decomp(csv_path: Path, out_file: Path):
    """Plots top contributing residues from gmx_MMPBSA decomp CSV."""
    if not csv_path.exists(): return
    try:
        res_data = []
        with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = [p.strip() for p in line.split(",")]
                # Decomp row format: R:Chain:ResName:ResID, ..., TOTAL
                if len(parts) >= 7 and parts[0].startswith(("R:", "L:")):
                    res_tag = parts[0]
                    if "UNK" in res_tag or "LIG" in res_tag: continue
                    val_str = parts[-1].split("+/-")[0].split("±")[0].strip()
                    try:
                        val = float(val_str)
                        if abs(val) >= 0.5:
                            p_parts = res_tag.split(":")
                            label = f"{p_parts[1]}:{p_parts[2]}{p_parts[3]}" if len(p_parts) >= 4 else res_tag
                            res_data.append({"Residue": label, "Energy": val})
                    except ValueError:
                        continue

        if not res_data: return
        df_plot = pd.DataFrame(res_data).sort_values("Energy")

        plt.figure(figsize=(10, max(4, len(df_plot) * 0.25)))
        colors = ["#16a34a" if x < 0 else "#dc2626" for x in df_plot["Energy"]]
        plt.barh(df_plot["Residue"], df_plot["Energy"], color=colors, edgecolor="black")
        plt.axvline(x=0, color='black', linewidth=1)
        plt.xlabel("ΔG Contribution (kcal/mol)")
        plt.title("Per-Residue MM-PBSA Energy Decomposition (Hot-Spots)")
        plt.grid(axis='x', alpha=0.3)
        plt.tight_layout()
        plt.savefig(out_file, dpi=300)
        plt.close()
    except Exception as e:
        import logging
        logging.warning(f"Failed to plot MM-PBSA decomp: {e}")