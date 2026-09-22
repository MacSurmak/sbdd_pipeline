"""
File: src/utils/md_metrics.py
Description: Clean mathematical extractors for MD metrics and statistical inefficiency.
"""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.lib.distances import capped_distance
from loguru import logger


def compute_rmsd_rmsf(pdb_path: str, xtc_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculates Protein Backbone RMSD, Ligand Heavy RMSD (drift), and Protein C-alpha RMSF."""
    u = mda.Universe(pdb_path, xtc_path)
    ref = mda.Universe(pdb_path)
    
    prot_bb = u.select_atoms("protein and backbone")
    ref_bb = ref.select_atoms("protein and backbone")
    lig_heavy = u.select_atoms("resname UNK and not name H*")
    ref_lig = ref.select_atoms("resname UNK and not name H*")
    
    rmsd_results = []
    for ts in u.trajectory:
        # 1. Align the entire system based on the protein backbone
        align.alignto(u.atoms, ref.atoms, select="protein and backbone")
        time_ns = ts.time / 1000.0
        
        # 2. Calculate RMSD. Ligand RMSD tracks physical unbinding (drift)
        prot_rmsd = rms.rmsd(prot_bb.positions, ref_bb.positions)
        lig_rmsd = rms.rmsd(lig_heavy.positions, ref_lig.positions) if lig_heavy.n_atoms > 0 else np.nan
        
        rmsd_results.append({"Time_ns": time_ns, "Prot_BB_RMSD": prot_rmsd, "Lig_Heavy_RMSD": lig_rmsd})
        
    df_rmsd = pd.DataFrame(rmsd_results)
    
    calphas = u.select_atoms("protein and name CA")
    R = rms.RMSF(calphas).run()
    
    df_rmsf = pd.DataFrame({
        "ResID": calphas.resids,
        "ResName": calphas.resnames,
        "RMSF": R.results.rmsf
    })
    
    return df_rmsd, df_rmsf

def compute_contacts(pdb_path: str, xtc_path: str, cutoffs: dict) -> pd.DataFrame:
    """Computes dynamic contacts occupancy (Total, Hydrophobic, Polar)."""
    u = mda.Universe(pdb_path, xtc_path)
    lig = u.select_atoms("resname UNK")
    prot = u.select_atoms("protein")
    
    if lig.n_atoms == 0 or prot.n_atoms == 0:
        return pd.DataFrame()
        
    l_heavy, p_heavy = lig.select_atoms("not name H*"), prot.select_atoms("not name H*")
    l_c, p_c = lig.select_atoms("name C*"), prot.select_atoms("name C*")
    l_pol = lig.select_atoms("name N* or name O* or name S* or name P*")
    p_pol = prot.select_atoms("name N* or name O* or name S* or name P*")
    
    unique_resids = np.unique(p_heavy.resids)
    counts = {r: {"Total": 0, "Hydrophobic": 0, "Polar": 0} for r in unique_resids}
    n_frames = len(u.trajectory)
    
    c_heavy, c_hyd, c_pol = cutoffs.get("heavy", 3.5), cutoffs.get("hydrophobic", 4.0), cutoffs.get("polar", 3.4)
    
    for _ in u.trajectory:
        # Total heavy contacts
        pairs_tot = capped_distance(l_heavy.positions, p_heavy.positions, max_cutoff=c_heavy, return_distances=False)
        if len(pairs_tot) > 0:
            for r in np.unique(p_heavy.resids[pairs_tot[:, 1]]): counts[r]["Total"] += 1
                
        # Hydrophobic (C-C)
        if l_c.n_atoms > 0 and p_c.n_atoms > 0:
            pairs_hyd = capped_distance(l_c.positions, p_c.positions, max_cutoff=c_hyd, return_distances=False)
            if len(pairs_hyd) > 0:
                for r in np.unique(p_c.resids[pairs_hyd[:, 1]]): counts[r]["Hydrophobic"] += 1
                    
        # Polar (N/O/S/P)
        if l_pol.n_atoms > 0 and p_pol.n_atoms > 0:
            pairs_pol = capped_distance(l_pol.positions, p_pol.positions, max_cutoff=c_pol, return_distances=False)
            if len(pairs_pol) > 0:
                for r in np.unique(p_pol.resids[pairs_pol[:, 1]]): counts[r]["Polar"] += 1

    results = []
    for r, c in counts.items():
        occ_tot = (c["Total"]/n_frames)*100.0
        if occ_tot > 1.0:
            results.append({
                "ResID": r, "ResName": prot.select_atoms(f"resid {r}").resnames[0],
                "Total_%": occ_tot, "Hydrophobic_%": (c["Hydrophobic"]/n_frames)*100.0,
                "Polar_%": (c["Polar"]/n_frames)*100.0
            })
            
    return pd.DataFrame(results).sort_values("Total_%", ascending=False)