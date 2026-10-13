#!/usr/bin/env python3
"""
File: bin/md_metrics.py
Description: Standalone CLI to calculate MD trajectory metrics (RMSD, RMSF, Contacts).
Outputs clean CSV tables and JSON summary without graphical dependencies.
"""
import argparse
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import capped_distance


def compute_rmsd_rmsf(
    struct_path: str, xtc_path: str, ligand_resname: str, cofactor_resnames: list
):
    # Load trajectory using GROMACS GRO or topology file natively without PDB column corruption
    u = mda.Universe(struct_path, xtc_path)
    ref = mda.Universe(struct_path)

    prot_bb = u.select_atoms("protein and backbone")
    ref_bb = ref.select_atoms("protein and backbone")
    prot_heavy = u.select_atoms("protein and not name H*")
    ref_prot_heavy = ref.select_atoms("protein and not name H*")

    # Strict declarative targets (Fail-Fast if residue not present)
    targets = []
    if ligand_resname:
        lig_sel = u.select_atoms(f"resname {ligand_resname} and not name H*")
        ref_lig = ref.select_atoms(f"resname {ligand_resname} and not name H*")
        if lig_sel.n_atoms == 0:
            lig_sel = u.select_atoms(f"resname {ligand_resname}")
            ref_lig = ref.select_atoms(f"resname {ligand_resname}")

        if lig_sel.n_atoms == 0:
            available = [
                r.resname
                for r in u.residues
                if r.atoms.select_atoms("protein").n_atoms == 0
            ]
            raise ValueError(
                f"FATAL: Ligand residue '{ligand_resname}' not found in structure! "
                f"Available non-protein residues: {set(available)}"
            )

        targets.append((f"Lig_{ligand_resname}_RMSD_A", lig_sel, ref_lig))

    for cof in cofactor_resnames:
        cof = cof.strip()
        if not cof:
            continue
        cof_sel = u.select_atoms(f"resname {cof} and not name H*")
        ref_cof = ref.select_atoms(f"resname {cof} and not name H*")
        if cof_sel.n_atoms == 0:
            # Monoatomic ions (MG, ZN, CA)
            cof_sel = u.select_atoms(f"resname {cof}")
            ref_cof = ref.select_atoms(f"resname {cof}")
        if cof_sel.n_atoms > 0:
            targets.append((f"Cofactor_{cof}_RMSD_A", cof_sel, ref_cof))

    rmsd_results = []
    for ts in u.trajectory:
        align.alignto(u.atoms, ref.atoms, select="protein and backbone")
        time_ns = ts.time / 1000.0

        entry = {
            "Time_ns": round(time_ns, 3),
            "Prot_BB_RMSD_A": round(
                float(rms.rmsd(prot_bb.positions, ref_bb.positions)), 3
            ),
            "Prot_Heavy_RMSD_A": round(
                float(rms.rmsd(prot_heavy.positions, ref_prot_heavy.positions)), 3
            ),
        }

        for col_name, u_sel, ref_sel in targets:
            r_val = float(rms.rmsd(u_sel.positions, ref_sel.positions))
            entry[col_name] = round(r_val, 3)

        rmsd_results.append(entry)

    df_rmsd = pd.DataFrame(rmsd_results)

    calphas = u.select_atoms("protein and name CA")
    R = rms.RMSF(calphas).run()

    df_rmsf = pd.DataFrame(
        {
            "ResID": calphas.resids,
            "ResName": calphas.resnames,
            "RMSF_A": [round(float(v), 3) for v in R.results.rmsf],
        }
    )

    return df_rmsd, df_rmsf


def compute_contacts(
    struct_path: str,
    xtc_path: str,
    ligand_resname: str,
    c_heavy: float,
    c_hyd: float,
    c_pol: float,
):
    u = mda.Universe(struct_path, xtc_path)
    prot = u.select_atoms("protein")

    if prot.n_atoms == 0 or not ligand_resname:
        return pd.DataFrame()

    lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
    if lig.n_atoms == 0:
        lig = u.select_atoms(f"resname {ligand_resname}")

    if lig.n_atoms == 0:
        return pd.DataFrame()
    l_heavy = lig.select_atoms("not name H*")
    p_heavy = prot.select_atoms("not name H*")
    l_c = lig.select_atoms("name C*")
    p_c = prot.select_atoms("name C*")
    l_pol = lig.select_atoms("name N* or name O* or name S* or name P*")
    p_pol = prot.select_atoms("name N* or name O* or name S* or name P*")

    unique_resids = np.unique(p_heavy.resids)
    counts = {r: {"Total": 0, "Hydrophobic": 0, "Polar": 0} for r in unique_resids}
    n_frames = len(u.trajectory)

    for _ in u.trajectory:
        pairs_tot = capped_distance(
            l_heavy.positions,
            p_heavy.positions,
            max_cutoff=c_heavy,
            return_distances=False,
        )
        if len(pairs_tot) > 0:
            for r in np.unique(p_heavy.resids[pairs_tot[:, 1]]):
                counts[r]["Total"] += 1

        if l_c.n_atoms > 0 and p_c.n_atoms > 0:
            pairs_hyd = capped_distance(
                l_c.positions, p_c.positions, max_cutoff=c_hyd, return_distances=False
            )
            if len(pairs_hyd) > 0:
                for r in np.unique(p_c.resids[pairs_hyd[:, 1]]):
                    counts[r]["Hydrophobic"] += 1

        if l_pol.n_atoms > 0 and p_pol.n_atoms > 0:
            pairs_pol = capped_distance(
                l_pol.positions,
                p_pol.positions,
                max_cutoff=c_pol,
                return_distances=False,
            )
            if len(pairs_pol) > 0:
                for r in np.unique(p_pol.resids[pairs_pol[:, 1]]):
                    counts[r]["Polar"] += 1

    results = []
    for r, c in counts.items():
        occ_tot = (c["Total"] / n_frames) * 100.0
        if occ_tot >= 1.0:
            res_atoms = prot.select_atoms(f"resid {r}")
            res_name = res_atoms.resnames[0] if len(res_atoms) > 0 else "UNK"
            results.append(
                {
                    "ResID": int(r),
                    "ResName": res_name,
                    "Total_Occ_Pct": round(occ_tot, 1),
                    "Hydrophobic_Occ_Pct": round(
                        (c["Hydrophobic"] / n_frames) * 100.0, 1
                    ),
                    "Polar_Occ_Pct": round((c["Polar"] / n_frames) * 100.0, 1),
                }
            )

    df_contacts = pd.DataFrame(results)
    if not df_contacts.empty:
        df_contacts = df_contacts.sort_values("Total_Occ_Pct", ascending=False)
    return df_contacts


def generate_rmsd_plot(df_rmsd: pd.DataFrame, out_png: Path):
    """Renders a clean multi-component RMSD plot without graphical display."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.figure(figsize=(10, 5), dpi=300)
        time_col = df_rmsd["Time_ns"]

        if "Prot_BB_RMSD_A" in df_rmsd.columns:
            plt.plot(
                time_col,
                df_rmsd["Prot_BB_RMSD_A"],
                label="Protein Backbone",
                color="#1f77b4",
                linewidth=1.5,
            )
        if "Prot_Heavy_RMSD_A" in df_rmsd.columns:
            plt.plot(
                time_col,
                df_rmsd["Prot_Heavy_RMSD_A"],
                label="Protein Heavy Atoms",
                color="#aec7e8",
                linewidth=1.0,
                alpha=0.7,
            )

        for col in df_rmsd.columns:
            if col.startswith("Lig_"):
                label_name = col.replace("Lig_", "").replace("_RMSD_A", " (Ligand)")
                plt.plot(
                    time_col,
                    df_rmsd[col],
                    label=label_name,
                    color="#d62728",
                    linewidth=2.0,
                )
            elif col.startswith("Cofactor_"):
                label_name = col.replace("Cofactor_", "").replace(
                    "_RMSD_A", " (Cofactor)"
                )
                plt.plot(
                    time_col,
                    df_rmsd[col],
                    label=label_name,
                    color="#ff7f0e",
                    linestyle="--",
                    linewidth=1.6,
                )

        # Reference threshold line for stability
        plt.axhline(
            y=2.0,
            color="#7f7f7f",
            linestyle=":",
            linewidth=1.0,
            alpha=0.7,
            label="2.0 Å Threshold",
        )

        plt.title("Molecular Dynamics RMSD Trajectory", fontsize=12, fontweight="bold")
        plt.xlabel("Simulation Time (ns)", fontsize=10)
        plt.ylabel("RMSD (Å)", fontsize=10)
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.legend(loc="upper left", frameon=True)
        plt.tight_layout()
        plt.savefig(out_png)
        plt.close()
    except Exception as e:
        print(f"Warning: Could not render RMSD plot: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate MD metrics (RMSD, RMSF, Contacts)"
    )
    parser.add_argument(
        "--struct", required=True, help="Clean reference structure in GROMACS PDB/GRO format"
    )
    parser.add_argument("--xtc", required=True, help="Clean Centered Trajectory XTC")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--ligand", default="", help="Canonical ligand residue name (e.g. KVR)"
    )
    parser.add_argument(
        "--cofactors",
        nargs="*",
        default=[],
        help="Cofactor residue names (e.g. ATP MG)",
    )
    parser.add_argument("--heavy_cutoff", type=float, required=True)
    parser.add_argument("--hydrophobic_cutoff", type=float, required=True)
    parser.add_argument("--polar_cutoff", type=float, required=True)
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        df_rmsd, df_rmsf = compute_rmsd_rmsf(
            args.struct, args.xtc, args.ligand, args.cofactors
        )
        df_rmsd.to_csv(out_dir / "rmsd.csv", index=False)
        df_rmsf.to_csv(out_dir / "rmsf.csv", index=False)
        generate_rmsd_plot(df_rmsd, out_dir / "rmsd_plot.png")

        df_contacts = compute_contacts(
            args.struct,
            args.xtc,
            args.ligand,
            c_heavy=args.heavy_cutoff,
            c_hyd=args.hydrophobic_cutoff,
            c_pol=args.polar_cutoff,
        )
        df_contacts.to_csv(out_dir / "contacts.csv", index=False)

        summary = {
            "Prot_BB_RMSD_Mean_A": round(float(df_rmsd["Prot_BB_RMSD_A"].mean()), 3),
            "Prot_BB_RMSD_Std_A": round(float(df_rmsd["Prot_BB_RMSD_A"].std()), 3),
            "Prot_BB_RMSD_Last_A": round(float(df_rmsd["Prot_BB_RMSD_A"].iloc[-1]), 3),
            "Molecules": {},
        }

        # Dynamically record stats for each ligand and cofactor
        first_lig = None
        for col in df_rmsd.columns:
            if col.startswith("Lig_") or col.startswith("Cofactor_"):
                if col.startswith("Lig_") and first_lig is None:
                    first_lig = col
                s = df_rmsd[col].dropna()
                if not s.empty:
                    summary["Molecules"][col] = {
                        "Mean_A": round(float(s.mean()), 3),
                        "Std_A": round(float(s.std()), 3),
                        "Last_A": round(float(s.iloc[-1]), 3),
                    }

        # Backwards-compatible fields
        if first_lig:
            summary["Lig_Heavy_RMSD_Mean_A"] = summary["Molecules"][first_lig]["Mean_A"]
            summary["Lig_Heavy_RMSD_Last_A"] = summary["Molecules"][first_lig]["Last_A"]
        else:
            summary["Lig_Heavy_RMSD_Mean_A"] = None
            summary["Lig_Heavy_RMSD_Last_A"] = None

        summary["Top_Contacts"] = (
            df_contacts.head(10).to_dict(orient="records")
            if not df_contacts.empty
            else []
        )

        with open(out_dir / "md_metrics_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        print(f"MD metrics calculated successfully. Saved to {out_dir}")

    except Exception as e:
        print(f"Error calculating MD metrics: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
