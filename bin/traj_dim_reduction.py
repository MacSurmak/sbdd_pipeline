#!/usr/bin/env python3
"""
bin/traj_dim_reduction.py
Atomic Stage 1: In-memory coordinate streaming, on-the-fly Kabsch superimposition,
and essential dynamics dimensionality reduction (PCA and TICA + 2D UMAP).
Completely avoids writing intermediate temporary XTC files to disk.
"""

import argparse
import json
import sys
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix

try:
    import umap
except ImportError:
    print("FATAL: umap-learn is required in current environment.", file=sys.stderr)
    sys.exit(1)


def parse_ndx(ndx_path: Path) -> dict[str, list[int]]:
    """Parses standard GROMACS NDX file into dictionary of group_name -> 0-based atom indices."""
    groups = {}
    current_group = None
    with open(ndx_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith("[") and line_str.endswith("]"):
                current_group = line_str[1:-1].strip()
                groups[current_group] = []
            elif current_group is not None:
                # GROMACS NDX atom indices are 1-based; convert to 0-based for MDAnalysis
                groups[current_group].extend([int(x) - 1 for x in line_str.split()])
    return groups


def resolve_atom_group(
    u: mda.Universe,
    group_name: str,
    ndx_groups: dict[str, list[int]],
    ligand_resname: str,
    pocket_radius: float,
    desc: str,
) -> mda.AtomGroup:
    """
    Deterministically resolves atom selection via NDX, macros, or explicit queries.
    Strict Fail-Fast law: crashes immediately if selection resolves to 0 atoms.
    """
    g_clean = group_name.strip()
    g_lower = g_clean.lower()

    # 1. Match named group from index.ndx
    if ndx_groups:
        for ndx_name, idx_list in ndx_groups.items():
            if ndx_name.strip().lower() == g_lower:
                ag = u.atoms[idx_list]
                if len(ag) > 0:
                    return ag

    # 2. Match standard biophysical macros
    if g_lower in ["backbone", "bb"]:
        ag = u.select_atoms("protein and (name CA or name C or name N)")
        if len(ag) > 0:
            return ag
    elif g_lower in ["c-alpha", "ca"]:
        ag = u.select_atoms("protein and name CA")
        if len(ag) > 0:
            return ag
    elif g_lower in ["pocket"]:
        if not ligand_resname:
            raise ValueError(
                f"FATAL: Group '{group_name}' for {desc} requires an explicit ligand residue name, "
                "but none was provided and system is Apo."
            )
        lig_atoms = u.select_atoms(f"resname {ligand_resname} and not name H*")
        if len(lig_atoms) == 0:
            raise ValueError(
                f"FATAL: Reference ligand '{ligand_resname}' not found in structure for {desc} pocket definition!"
            )
        ag = u.select_atoms(
            f"(same residue as (protein and not name H* and not type H and around {pocket_radius} resname {ligand_resname})) and not name H* and not type H"
        )
        if len(ag) > 0:
            return ag
    elif g_lower in ["ligand", "lig"]:
        if not ligand_resname:
            raise ValueError(
                f"FATAL: Group '{group_name}' for {desc} requires a ligand residue name, but none was provided."
            )
        ag = u.select_atoms(f"resname {ligand_resname} and not name H*")
        if len(ag) > 0:
            return ag

    # 3. Match raw MDAnalysis selection syntax
    try:
        ag = u.select_atoms(g_clean)
        if len(ag) > 0:
            return ag
    except Exception as e:
        raise ValueError(
            f"FATAL: Failed to evaluate selection string '{g_clean}' for {desc}: {e}"
        )

    avail_ndx = list(ndx_groups.keys()) if ndx_groups else []
    raise ValueError(
        f"FATAL: Selection group '{group_name}' for {desc} resolved to 0 atoms!\n"
        f"Available groups in index.ndx: {avail_ndx}"
    )


def compute_tica(
    rep_coords_list: list[np.ndarray], lag_steps: int, variance_cutoff: float
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Computes Time-lagged Independent Component Analysis (TICA) across multi-replicate ensembles.
    Enforces strict physical boundaries: zero cross-correlation across independent trajectories.
    Eliminates rank deficiency via SVD whitening before solving standard Hermitian eigensystem.
    """
    from scipy import linalg

    total_samples = sum(len(r) for r in rep_coords_list)
    n_features = rep_coords_list[0].shape[1]

    # Global centroid across all replicates
    all_coords = np.vstack(rep_coords_list)
    mean_vec = np.mean(all_coords, axis=0)
    all_centered = all_coords - mean_vec

    # 1. PCA whitening to eliminate rank deficiency and collinearity in Cartesian coordinates
    # X = U * S * Vt -> Covariance = (1 / (N - 1)) * V * S^2 * Vt
    u_mat, s_vals, vt_mat = linalg.svd(all_centered, full_matrices=False)
    var_all = (s_vals ** 2) / float(max(1, total_samples - 1))
    sum_var = np.sum(var_all)
    cum_var = np.cumsum(var_all) / sum_var if sum_var > 0 else np.zeros_like(var_all)

    # Retain components covering variance_cutoff without singular modes (singular value > 1e-4)
    valid_sing = s_vals > 1e-4
    m_cand = int(np.searchsorted(cum_var, variance_cutoff) + 1)
    m_components = max(2, min(m_cand, int(np.sum(valid_sing)), 50, total_samples - 1))

    v_white = vt_mat[:m_components, :].T  # Shape: (n_features, m_components)
    s_white = s_vals[:m_components]       # Singular values for scaling

    # Whitening operator: transforms coordinates such that C(0) = Identity_M
    # y = x * V * diag(sqrt(N - 1) / S)
    white_scale = np.sqrt(float(max(1, total_samples - 1))) / s_white
    w_operator = v_white * white_scale[np.newaxis, :]

    # 2. Accumulate C(tau) strictly within each individual replicate trajectory
    c_tau_accum = np.zeros((m_components, m_components), dtype=np.float64)
    total_tau_pairs = 0

    whitened_reps = []
    for r_coords in rep_coords_list:
        if len(r_coords) <= lag_steps:
            continue
        r_centered = r_coords - mean_vec
        r_whitened = np.dot(r_centered, w_operator)
        whitened_reps.append(r_whitened)

        y_t0 = r_whitened[:-lag_steps]
        y_tau = r_whitened[lag_steps:]
        c_tau_accum += np.dot(y_t0.T, y_tau)
        total_tau_pairs += len(y_t0)

    if total_tau_pairs == 0:
        raise ValueError(
            f"FATAL: All individual replicates are shorter than TICA lag steps ({lag_steps})! "
            "Cannot calculate time-lagged cross-correlations."
        )

    c_tau = c_tau_accum / float(total_tau_pairs)
    c_tau = 0.5 * (c_tau + c_tau.T)  # Enforce exact physical microscopic reversibility

    # Since C(0) is strictly Identity_M in whitened space, the generalized problem
    # C(tau) * v = lambda * C(0) * v reduces to a perfectly conditioned standard eigensystem:
    eigvals, eigvecs = linalg.eigh(c_tau)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    pos_vals = np.maximum(eigvals, 0.0)
    sum_vals = np.sum(pos_vals)
    var_exp = (pos_vals / sum_vals) if sum_vals > 0 else pos_vals

    k = max(2, min(10, m_components))
    all_whitened = np.dot(all_centered, w_operator)
    proj = np.dot(all_whitened, eigvecs[:, :k])

    return proj, var_exp, k


def compute_pca(
    data: np.ndarray, variance_cutoff: float
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Computes Principal Component Analysis (PCA) on Cartesian covariance matrix.
    Utilizes SVD decomposition on centered coordinates to eliminate memory blowups.
    """
    from scipy import linalg

    n_samples, n_features = data.shape
    mean_vec = np.mean(data, axis=0)
    centered = data - mean_vec

    _, s_vals, vt_mat = linalg.svd(centered, full_matrices=False)
    var_exp = (s_vals ** 2) / float(max(1, n_samples - 1))
    sum_var = np.sum(var_exp)
    if sum_var > 0:
        var_exp = var_exp / sum_var

    cum_var = np.cumsum(var_exp)
    k_cand = int(np.searchsorted(cum_var, variance_cutoff) + 1)
    k = max(2, min(10, k_cand, n_features, n_samples - 1))

    eigvecs = vt_mat.T
    proj = np.dot(centered, eigvecs[:, :k])
    return proj, var_exp, k


def main():
    parser = argparse.ArgumentParser(
        description="Stream trajectories, superimpose in memory, and reduce dimensionality"
    )
    parser.add_argument("--struct", required=True, help="Clean reference structure PDB")
    parser.add_argument("--xtc", nargs="+", required=True, help="Clean trajectory XTC file(s)")
    parser.add_argument("--rep_names", nargs="+", required=True, help="Explicit physical replicate names matching each XTC")
    parser.add_argument("--ndx", default="", help="Optional GROMACS index.ndx file")
    parser.add_argument("--align_group", required=True, help="Basis group for superimposition")
    parser.add_argument("--cluster_group", required=True, help="Feature group for coordinates")
    parser.add_argument("--ligand_resname", default="", help="Ligand residue name (if applicable)")
    parser.add_argument("--pocket_radius", type=float, required=True, help="Pocket radius (Å)")
    parser.add_argument("--dim_reduction", choices=["pca", "tica"], required=True)
    parser.add_argument("--tica_lag_ns", type=float, required=True, help="Lag time for TICA (ns)")
    parser.add_argument("--variance_cutoff", type=float, required=True, help="Variance cutoff fraction")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    struct_p = Path(args.struct).resolve()
    if not struct_p.exists() or struct_p.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Reference structure missing or empty: {struct_p}")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ndx_groups = parse_ndx(Path(args.ndx).resolve()) if args.ndx and Path(args.ndx).exists() else {}

    ref_u = mda.Universe(str(struct_p))
    ref_align_atoms = resolve_atom_group(
        ref_u, args.align_group, ndx_groups, args.ligand_resname, args.pocket_radius, "Alignment"
    )
    ref_cluster_atoms = resolve_atom_group(
        ref_u, args.cluster_group, ndx_groups, args.ligand_resname, args.pocket_radius, "Clustering Feature"
    )

    ref_align_pos = ref_align_atoms.positions.copy()
    ref_align_com = ref_align_pos.mean(axis=0)
    ref_align_centered = ref_align_pos - ref_align_com

    align_indices = ref_align_atoms.indices
    cluster_indices = ref_cluster_atoms.indices

    print(
        f"[DIM_REDUCTION] Reference loaded: {struct_p.name}\n"
        f"  - Alignment Group ('{args.align_group}'): {len(align_indices)} atoms\n"
        f"  - Feature Group ('{args.cluster_group}'): {len(cluster_indices)} atoms",
        flush=True,
    )

    clean_xtcs = []
    for raw_arg in args.xtc:
        for chunk in str(raw_arg).split():
            c_clean = chunk.strip()
            if c_clean:
                p_xtc = Path(c_clean).resolve()
                if not p_xtc.exists() or p_xtc.stat().st_size == 0:
                    raise FileNotFoundError(f"FATAL: Trajectory file missing or empty: {p_xtc}")
                clean_xtcs.append(p_xtc)

    clean_reps = []
    for raw_rep in args.rep_names:
        for r_chunk in str(raw_rep).split():
            r_c = r_chunk.strip()
            if r_c:
                clean_reps.append(r_c)

    if len(clean_xtcs) != len(clean_reps):
        raise ValueError(
            f"FATAL: Trajectory files count ({len(clean_xtcs)}) does not match explicit replicate names count ({len(clean_reps)})!"
        )

    coords_by_rep = []
    frame_time_map = []
    global_frame = 0
    dt_ps = None

    for xtc_p, rep_name in zip(clean_xtcs, clean_reps):
        u = mda.Universe(str(struct_p), str(xtc_p))

        if len(u.trajectory) < 2:
            raise ValueError(
                f"FATAL: Trajectory {xtc_p.name} for replicate '{rep_name}' has fewer than 2 frames! File is corrupted."
            )

        if dt_ps is None and len(u.trajectory) >= 2:
            dt_ps = float(u.trajectory[1].time - u.trajectory[0].time)

        rep_frames_coords = []
        for ts in u.trajectory:
            mob_align = u.atoms[align_indices].positions
            mob_align_com = mob_align.mean(axis=0)
            mob_align_centered = mob_align - mob_align_com

            r_mat, _ = rotation_matrix(mob_align_centered, ref_align_centered)
            mob_cluster = u.atoms[cluster_indices].positions
            aligned_cluster = np.dot(mob_cluster - mob_align_com, r_mat.T) + ref_align_com

            flat_coords = aligned_cluster.reshape(-1)
            rep_frames_coords.append(flat_coords)
            frame_time_map.append({
                "global_frame": global_frame,
                "replicate": rep_name,
                "time_ns": round(float(ts.time) / 1000.0, 4),
            })
            global_frame += 1

        coords_by_rep.append(np.array(rep_frames_coords, dtype=np.float64))

    all_coords_matrix = np.vstack(coords_by_rep)
    n_frames, n_features = all_coords_matrix.shape
    print(f"[DIM_REDUCTION] Assembled matrix: {n_frames} frames × {n_features} features across {len(clean_reps)} replicates.", flush=True)

    dt_ns = (dt_ps / 1000.0) if dt_ps and dt_ps > 0 else 0.01
    lag_steps = max(1, int(round(args.tica_lag_ns / dt_ns)))

    if args.dim_reduction == "pca":
        proj, var_exp, k = compute_pca(all_coords_matrix, args.variance_cutoff)
    else:
        proj, var_exp, k = compute_tica(coords_by_rep, lag_steps, args.variance_cutoff)

    print(
        f"[DIM_REDUCTION] {args.dim_reduction.upper()} completed: {k} components selected "
        f"(Variance: C1={var_exp[0]*100:.1f}%, C2={var_exp[1]*100:.1f}%).",
        flush=True,
    )

    n_neighbors = min(30, max(5, n_frames - 1))
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=0.1,
        metric="euclidean",
        n_jobs=1,
        random_state=42,
    )
    umap_2d = reducer.fit_transform(proj)

    np.save(outdir / "projections.npy", proj)
    np.save(outdir / "umap_2d.npy", umap_2d)
    np.save(outdir / "eigenvalues.npy", var_exp)
    np.save(outdir / "aligned_coords.npy", all_coords_matrix)

    with open(outdir / "frame_time_map.json", "w", encoding="utf-8") as f:
        json.dump(frame_time_map, f, indent=2)

    meta = {
        "align_group": args.align_group,
        "cluster_group": args.cluster_group,
        "dim_reduction": args.dim_reduction,
        "n_frames": n_frames,
        "n_features": n_features,
        "n_align_atoms": len(align_indices),
        "n_cluster_atoms": len(cluster_indices),
        "k_components": k,
        "variance_explained": [round(float(v), 4) for v in var_exp[:k]],
    }
    with open(outdir / "feature_meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(f"✓ Saved dimensionality reduction deliverables in {outdir}", flush=True)


if __name__ == "__main__":
    main()
