"""
File: cli/core/clustering.py
Description: Single Source of Truth (SSOT) for conformational ensemble clustering of docked poses
and trajectory sub-ensembles. Supports K >= 1 distance-based hierarchical clustering,
symmetry-aware vectorized RMSD matrices, and minimum-energy basin leader selection.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Final
import numpy as np
from rdkit import Chem
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

DEFAULT_POSE_CLUSTER_CUTOFF_A: Final[float] = 1.8
MAX_ENSEMBLE_MODES_PER_ENGINE: Final[int] = 16


@dataclass(frozen=True)
class PoseRecord:
    """Compact in-memory representation of a single docked pose."""

    engine: str
    mode_rank: int
    energy: float
    cnn_score: float | None
    cnn_affinity: float | None
    coords: np.ndarray


@dataclass(frozen=True)
class SubstanceClusterSummary:
    """Result of conformational clustering for a single substance in a receptor frame."""

    title: str
    n_glide_modes: int
    n_unidock_modes: int
    n_clusters: int
    winning_cluster_size: int
    is_consensus_basin: bool
    glide_best_mode: int | None
    glide_score: float | None
    glide_cnn_score: float | None
    glide_cnn_affinity: float | None
    unidock_best_mode: int | None
    unidock_vina: float | None
    unidock_cnn_score: float | None
    unidock_cnn_affinity: float | None
    pose_rmsd: float | None


def _build_untyped_heavy_graph(mol: Chem.Mol) -> Chem.Mol:
    """Builds a bond-order and charge-agnostic heavy-atom graph for symmetry perception."""
    m_heavy = Chem.RemoveHs(mol, sanitize=False)
    q_mol = Chem.Mol(m_heavy)
    for b in q_mol.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in q_mol.GetAtoms():
        a.SetFormalCharge(0)
        a.SetIsotope(0)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(False)
        a.SetNumRadicalElectrons(0)
    q_mol.UpdatePropertyCache(strict=False)
    return q_mol


def _extract_heavy_coords(mol: Chem.Mol) -> np.ndarray:
    """Extracts heavy-atom 3D coordinates into a compact (N_heavy, 3) float64 array."""
    m_heavy = Chem.RemoveHs(mol, sanitize=False)
    n_heavy = m_heavy.GetNumAtoms()
    conf = m_heavy.GetConformer()
    arr = np.empty((n_heavy, 3), dtype=np.float64)
    for idx in range(n_heavy):
        pos = conf.GetAtomPosition(idx)
        arr[idx, 0] = pos.x
        arr[idx, 1] = pos.y
        arr[idx, 2] = pos.z
    return arr


def _same_heavy_atom_order(q_ref: Chem.Mol, q_cand: Chem.Mol) -> bool:
    """Fast O(N) check whether two heavy-atom graphs share identical atom and bond indexing."""
    if q_ref.GetNumAtoms() != q_cand.GetNumAtoms():
        return False
    if q_ref.GetNumBonds() != q_cand.GetNumBonds():
        return False
    for idx in range(q_ref.GetNumAtoms()):
        if (
            q_ref.GetAtomWithIdx(idx).GetAtomicNum()
            != q_cand.GetAtomWithIdx(idx).GetAtomicNum()
        ):
            return False
    return True


def compute_all_pairs_symmetry_rmsd(
    coords_stack: np.ndarray,
    automorphisms: list[tuple[int, ...]],
) -> np.ndarray:
    """
    Computes an all-against-all (M, M) symmetry-aware heavy-atom RMSD matrix
    using vectorized NumPy broadcasting across all graph automorphisms.
    """
    n_poses = coords_stack.shape[0]
    n_heavy = coords_stack.shape[1]
    if n_poses <= 1 or n_heavy == 0:
        return np.zeros((n_poses, n_poses), dtype=np.float64)

    best_sq_dist = np.full((n_poses, n_poses), np.inf, dtype=np.float64)
    for perm in automorphisms:
        perm_coords = coords_stack[:, perm, :]
        diff = coords_stack[:, np.newaxis, :, :] - perm_coords[np.newaxis, :, :, :]
        sq_dist = np.sum(diff**2, axis=(2, 3)) / float(n_heavy)
        best_sq_dist = np.minimum(best_sq_dist, sq_dist)

    rmsd_mat = np.sqrt(np.maximum(best_sq_dist, 0.0))
    np.fill_diagonal(rmsd_mat, 0.0)
    return 0.5 * (rmsd_mat + rmsd_mat.T)


def partition_distance_matrix(
    rmsd_mat: np.ndarray,
    cutoff_angstrom: float,
) -> np.ndarray:
    """
    Partitions an (M, M) RMSD distance matrix into conformational basins (1..K)
    using average-linkage hierarchical clustering. Naturally returns K=1 when all
    poses lie within cutoff_angstrom.
    """
    n_poses = rmsd_mat.shape[0]
    if n_poses <= 1:
        return np.ones(n_poses, dtype=np.int32)
    condensed = squareform(rmsd_mat, checks=False)
    z_mat = linkage(condensed, method="average")
    labels = fcluster(z_mat, t=cutoff_angstrom, criterion="distance")
    return labels.astype(np.int32)


def _read_prop_float(mol: Chem.Mol, candidates: list[str]) -> float | None:
    prop_map = {p.strip().lower(): p for p in mol.GetPropNames()}
    for cand in candidates:
        c_low = cand.lower()
        if c_low in prop_map:
            raw = mol.GetProp(prop_map[c_low]).strip()
            if raw and raw.lower() != "nan":
                return float(raw.split()[0])
    return None


def index_engine_poses_from_sdf(
    sdf_path: Path,
    engine_name: str,
    fallback_scores: dict[str, list[float]],
    ref_graphs: dict[str, tuple[Chem.Mol, list[tuple[int, ...]]]],
) -> dict[str, list[PoseRecord]]:
    """
    Streams an engine's pose SDF file, aligns heavy-atom coordinates to a canonical
    per-substance reference graph, sorts modes by native energy, and retains up to
    MAX_ENSEMBLE_MODES_PER_ENGINE best poses per substance.
    """
    raw_by_title: dict[str, list[PoseRecord]] = {}
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        return raw_by_title

    mode_counters: dict[str, int] = {}
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    for rec_idx, mol in enumerate(suppl, start=1):
        if mol is None:
            raise ValueError(
                f"FATAL: Corrupted molecule record at index {rec_idx} in {sdf_path}!"
            )
        title = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        if not title and mol.HasProp("s_m_title"):
            title = mol.GetProp("s_m_title").strip()
        if not title:
            continue

        prev_cnt = mode_counters[title] if title in mode_counters else 0
        m_idx = prev_cnt + 1
        mode_counters[title] = m_idx

        if engine_name == "glide":
            energy = _read_prop_float(
                mol,
                ["r_i_docking_score", "docking_score", "r_i_glide_gscore", "gscore"],
            )
            if energy is None and title in fallback_scores:
                sc_list = fallback_scores[title]
                if m_idx <= len(sc_list):
                    energy = sc_list[m_idx - 1]
        else:
            energy = _read_prop_float(
                mol, ["unidock_vina_score", "minimizedAffinity", "affinity"]
            )

        if energy is None:
            continue

        cnn_sc = _read_prop_float(mol, ["CNNscore"])
        cnn_aff = _read_prop_float(mol, ["CNNaffinity"])

        q_cur = _build_untyped_heavy_graph(mol)
        coords = _extract_heavy_coords(mol)

        if title not in ref_graphs:
            matches = [tuple(m) for m in q_cur.GetSubstructMatches(q_cur, uniquify=False)]
            if not matches:
                matches = [tuple(range(q_cur.GetNumAtoms()))]
            ref_graphs[title] = (q_cur, matches)
        else:
            q_ref, _ = ref_graphs[title]
            if q_cur.GetNumAtoms() != q_ref.GetNumAtoms():
                continue
            if not _same_heavy_atom_order(q_ref, q_cur):
                match = q_cur.GetSubstructMatch(q_ref)
                if not match:
                    continue
                coords = coords[list(match), :]

        if title not in raw_by_title:
            raw_by_title[title] = []
        raw_by_title[title].append(
            PoseRecord(
                engine=engine_name,
                mode_rank=m_idx,
                energy=energy,
                cnn_score=cnn_sc,
                cnn_affinity=cnn_aff,
                coords=coords,
            )
        )

    trimmed: dict[str, list[PoseRecord]] = {}
    for title, p_list in raw_by_title.items():
        p_list.sort(key=lambda p: p.energy)
        kept = p_list[:MAX_ENSEMBLE_MODES_PER_ENGINE]
        trimmed[title] = [
            PoseRecord(
                engine=p_rec.engine,
                mode_rank=rank_i,
                energy=p_rec.energy,
                cnn_score=p_rec.cnn_score,
                cnn_affinity=p_rec.cnn_affinity,
                coords=p_rec.coords,
            )
            for rank_i, p_rec in enumerate(kept, start=1)
        ]

    return trimmed


def evaluate_substance_pose_ensemble(
    title: str,
    glide_poses: list[PoseRecord],
    unidock_poses: list[PoseRecord],
    automorphisms: list[tuple[int, ...]],
    cutoff_angstrom: float,
) -> SubstanceClusterSummary:
    """
    Clusters the combined (Glide + Uni-Dock) or single-engine pose ensemble for a substance,
    identifies consensus basins (or single-engine basins), and selects the minimum-energy
    leader pose(s) within the optimal basin.
    """
    n_g = len(glide_poses)
    n_u = len(unidock_poses)
    all_poses = glide_poses + unidock_poses
    n_total = len(all_poses)

    if n_total == 0:
        return SubstanceClusterSummary(
            title=title,
            n_glide_modes=0,
            n_unidock_modes=0,
            n_clusters=0,
            winning_cluster_size=0,
            is_consensus_basin=False,
            glide_best_mode=None,
            glide_score=None,
            glide_cnn_score=None,
            glide_cnn_affinity=None,
            unidock_best_mode=None,
            unidock_vina=None,
            unidock_cnn_score=None,
            unidock_cnn_affinity=None,
            pose_rmsd=None,
        )

    coords_stack = np.stack([p.coords for p in all_poses], axis=0)
    rmsd_mat = compute_all_pairs_symmetry_rmsd(coords_stack, automorphisms)
    labels = partition_distance_matrix(rmsd_mat, cutoff_angstrom)
    unique_cids = sorted(list(set(int(x) for x in labels)))
    n_clusters = len(unique_cids)

    if n_g > 0 and n_u > 0:
        consensus_candidates: list[
            tuple[float, int, PoseRecord, PoseRecord, float]
        ] = []
        for cid in unique_cids:
            c_indices = [i for i in range(n_total) if int(labels[i]) == cid]
            g_in_c = [i for i in c_indices if all_poses[i].engine == "glide"]
            u_in_c = [i for i in c_indices if all_poses[i].engine == "unidock"]
            if g_in_c and u_in_c:
                best_g_idx = min(g_in_c, key=lambda i: all_poses[i].energy)
                best_u_idx = min(u_in_c, key=lambda i: all_poses[i].energy)
                sub_mat = rmsd_mat[np.ix_(g_in_c, u_in_c)]
                min_pair_rmsd = float(np.min(sub_mat))
                comb_energy = all_poses[best_g_idx].energy + all_poses[best_u_idx].energy
                consensus_candidates.append(
                    (
                        comb_energy,
                        len(c_indices),
                        all_poses[best_g_idx],
                        all_poses[best_u_idx],
                        min_pair_rmsd,
                    )
                )

        if consensus_candidates:
            consensus_candidates.sort(key=lambda x: (x[0], x[4]))
            _, win_size, win_g, win_u, win_rmsd = consensus_candidates[0]
            return SubstanceClusterSummary(
                title=title,
                n_glide_modes=n_g,
                n_unidock_modes=n_u,
                n_clusters=n_clusters,
                winning_cluster_size=win_size,
                is_consensus_basin=True,
                glide_best_mode=win_g.mode_rank,
                glide_score=round(win_g.energy, 3),
                glide_cnn_score=win_g.cnn_score,
                glide_cnn_affinity=win_g.cnn_affinity,
                unidock_best_mode=win_u.mode_rank,
                unidock_vina=round(win_u.energy, 3),
                unidock_cnn_score=win_u.cnn_score,
                unidock_cnn_affinity=win_u.cnn_affinity,
                pose_rmsd=round(win_rmsd, 3),
            )

        cross_mat = rmsd_mat[:n_g, n_g:]
        min_flat_idx = int(np.argmin(cross_mat))
        g_pair_i, u_pair_j = divmod(min_flat_idx, n_u)
        closest_rmsd = float(cross_mat[g_pair_i, u_pair_j])

        win_g = glide_poses[0]
        win_u = unidock_poses[0]
        win_g_cid = int(labels[0])
        win_size = sum(1 for x in labels if int(x) == win_g_cid)
        return SubstanceClusterSummary(
            title=title,
            n_glide_modes=n_g,
            n_unidock_modes=n_u,
            n_clusters=n_clusters,
            winning_cluster_size=win_size,
            is_consensus_basin=False,
            glide_best_mode=win_g.mode_rank,
            glide_score=round(win_g.energy, 3),
            glide_cnn_score=win_g.cnn_score,
            glide_cnn_affinity=win_g.cnn_affinity,
            unidock_best_mode=win_u.mode_rank,
            unidock_vina=round(win_u.energy, 3),
            unidock_cnn_score=win_u.cnn_score,
            unidock_cnn_affinity=win_u.cnn_affinity,
            pose_rmsd=round(closest_rmsd, 3),
        )

    e_min = min(p.energy for p in all_poses)
    basin_candidates: list[tuple[float, int, PoseRecord]] = []
    for cid in unique_cids:
        c_indices = [i for i in range(n_total) if int(labels[i]) == cid]
        best_idx = min(c_indices, key=lambda i: all_poses[i].energy)
        leader = all_poses[best_idx]
        if leader.energy <= e_min + 1.5:
            cnn_bonus = (
                2.0 * leader.cnn_score if leader.cnn_score is not None else 0.0
            )
            pop_bonus = 0.08 * float(len(c_indices))
            quality = -leader.energy + cnn_bonus + pop_bonus
            basin_candidates.append((quality, len(c_indices), leader))

    basin_candidates.sort(key=lambda x: x[0], reverse=True)
    _, win_size, win_pose = basin_candidates[0]

    if n_g > 0:
        return SubstanceClusterSummary(
            title=title,
            n_glide_modes=n_g,
            n_unidock_modes=0,
            n_clusters=n_clusters,
            winning_cluster_size=win_size,
            is_consensus_basin=False,
            glide_best_mode=win_pose.mode_rank,
            glide_score=round(win_pose.energy, 3),
            glide_cnn_score=win_pose.cnn_score,
            glide_cnn_affinity=win_pose.cnn_affinity,
            unidock_best_mode=None,
            unidock_vina=None,
            unidock_cnn_score=None,
            unidock_cnn_affinity=None,
            pose_rmsd=None,
        )

    return SubstanceClusterSummary(
        title=title,
        n_glide_modes=0,
        n_unidock_modes=n_u,
        n_clusters=n_clusters,
        winning_cluster_size=win_size,
        is_consensus_basin=False,
        glide_best_mode=None,
        glide_score=None,
        glide_cnn_score=None,
        glide_cnn_affinity=None,
        unidock_best_mode=win_pose.mode_rank,
        unidock_vina=round(win_pose.energy, 3),
        unidock_cnn_score=win_pose.cnn_score,
        unidock_cnn_affinity=win_pose.cnn_affinity,
        pose_rmsd=None,
    )
