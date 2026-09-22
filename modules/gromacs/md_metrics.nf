/*
 * modules/gromacs/md_metrics.nf
 * Calculates trajectory geometry metrics (RMSD, RMSF, Contacts)
 * Self-contained worker process: zero external script file staging.
 */

process MD_METRICS {
    tag "${meta.id}_rep${meta.rep}"
    label "process_medium"

    input:
    tuple val(meta), path(clean_xtc), path(ref_pdb)

    output:
    tuple val(meta), path("rmsd.csv"), path("rmsf.csv"), path("contacts.csv"), emit: tables
    tuple val(meta), path("md_metrics_summary.json"), emit: summary_json

    script:
    """
    cat << 'EOF' > md_metrics.py
import argparse, sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import capped_distance

def compute_rmsd_rmsf(pdb_path: str, xtc_path: str):
    u = mda.Universe(pdb_path, xtc_path)
    ref = mda.Universe(pdb_path)
    prot_bb = u.select_atoms("protein and backbone")
    ref_bb = ref.select_atoms("protein and backbone")
    lig_heavy = u.select_atoms("resname UNK and not name H*")
    ref_lig = ref.select_atoms("resname UNK and not name H*")
    
    rmsd_results = []
    for ts in u.trajectory:
        align.alignto(u.atoms, ref.atoms, select="protein and backbone")
        time_ns = ts.time / 1000.0
        prot_rmsd = float(rms.rmsd(prot_bb.positions, ref_bb.positions))
        lig_rmsd = float(rms.rmsd(lig_heavy.positions, ref_lig.positions)) if lig_heavy.n_atoms > 0 else np.nan
        rmsd_results.append({
            "Time_ns": round(time_ns, 3),
            "Prot_BB_RMSD_A": round(prot_rmsd, 3),
            "Lig_Heavy_RMSD_A": round(lig_rmsd, 3) if not np.isnan(lig_rmsd) else None
        })
        
    df_rmsd = pd.DataFrame(rmsd_results)
    calphas = u.select_atoms("protein and name CA")
    R = rms.RMSF(calphas).run()
    df_rmsf = pd.DataFrame({
        "ResID": calphas.resids,
        "ResName": calphas.resnames,
        "RMSF_A": [round(float(v), 3) for v in R.results.rmsf]
    })
    return df_rmsd, df_rmsf

def compute_contacts(pdb_path: str, xtc_path: str, c_heavy=3.5, c_hyd=4.0, c_pol=3.4):
    u = mda.Universe(pdb_path, xtc_path)
    lig = u.select_atoms("resname UNK")
    prot = u.select_atoms("protein")
    if lig.n_atoms == 0 or prot.n_atoms == 0:
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
        pairs_tot = capped_distance(l_heavy.positions, p_heavy.positions, max_cutoff=c_heavy, return_distances=False)
        if len(pairs_tot) > 0:
            for r in np.unique(p_heavy.resids[pairs_tot[:, 1]]):
                counts[r]["Total"] += 1
        if l_c.n_atoms > 0 and p_c.n_atoms > 0:
            pairs_hyd = capped_distance(l_c.positions, p_c.positions, max_cutoff=c_hyd, return_distances=False)
            if len(pairs_hyd) > 0:
                for r in np.unique(p_c.resids[pairs_hyd[:, 1]]):
                    counts[r]["Hydrophobic"] += 1
        if l_pol.n_atoms > 0 and p_pol.n_atoms > 0:
            pairs_pol = capped_distance(l_pol.positions, p_pol.positions, max_cutoff=c_pol, return_distances=False)
            if len(pairs_pol) > 0:
                for r in np.unique(p_pol.resids[pairs_pol[:, 1]]):
                    counts[r]["Polar"] += 1

    results = []
    for r, c in counts.items():
        occ_tot = (c["Total"] / n_frames) * 100.0
        if occ_tot >= 1.0:
            res_atoms = prot.select_atoms(f"resid {r}")
            res_name = res_atoms.resnames[0] if len(res_atoms) > 0 else "UNK"
            results.append({
                "ResID": int(r),
                "ResName": res_name,
                "Total_Occ_Pct": round(occ_tot, 1),
                "Hydrophobic_Occ_Pct": round((c["Hydrophobic"] / n_frames) * 100.0, 1),
                "Polar_Occ_Pct": round((c["Polar"] / n_frames) * 100.0, 1)
            })
    df_contacts = pd.DataFrame(results)
    if not df_contacts.empty:
        df_contacts = df_contacts.sort_values("Total_Occ_Pct", ascending=False)
    return df_contacts

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", required=True)
parser.add_argument("--xtc", required=True)
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
out_dir = Path(args.outdir)
out_dir.mkdir(parents=True, exist_ok=True)

df_rmsd, df_rmsf = compute_rmsd_rmsf(args.pdb, args.xtc)
df_rmsd.to_csv(out_dir / "rmsd.csv", index=False)
df_rmsf.to_csv(out_dir / "rmsf.csv", index=False)
df_contacts = compute_contacts(args.pdb, args.xtc)
df_contacts.to_csv(out_dir / "contacts.csv", index=False)
summary = {
    "Prot_BB_RMSD_Mean_A": round(float(df_rmsd["Prot_BB_RMSD_A"].mean()), 3),
    "Prot_BB_RMSD_Std_A": round(float(df_rmsd["Prot_BB_RMSD_A"].std()), 3),
    "Prot_BB_RMSD_Last_A": round(float(df_rmsd["Prot_BB_RMSD_A"].iloc[-1]), 3),
    "Lig_Heavy_RMSD_Mean_A": round(float(df_rmsd["Lig_Heavy_RMSD_A"].dropna().mean()), 3) if not df_rmsd["Lig_Heavy_RMSD_A"].dropna().empty else None,
    "Lig_Heavy_RMSD_Last_A": round(float(df_rmsd["Lig_Heavy_RMSD_A"].dropna().iloc[-1]), 3) if not df_rmsd["Lig_Heavy_RMSD_A"].dropna().empty else None,
    "Top_Contacts": df_contacts.head(10).to_dict(orient="records") if not df_contacts.empty else []
}
with open(out_dir / "md_metrics_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
EOF

    python3 md_metrics.py \
        --pdb ${ref_pdb} \
        --xtc ${clean_xtc} \
        --outdir .
    """

    stub:
    """
    touch rmsd.csv
    touch rmsf.csv
    touch contacts.csv
    touch md_metrics_summary.json
    """
}
