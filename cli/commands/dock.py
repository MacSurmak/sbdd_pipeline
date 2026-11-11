"""
File: cli/commands/dock.py
Description: CLI handler for 'sbdd dock'. Reads DB to construct Cross-Docking plan.
"""

import csv
import json
import math
import statistics
from pathlib import Path
import sys
from rdkit import Chem
from rich.console import Console
from rich.table import Table
from sqlalchemy.orm import Session
from cli.core.chem import (
    compute_pairwise_symmetry_rmsd_matrix,
    compute_symmetry_rmsd,
    get_microstate_identity,
)
from cli.core.config_loader import build_nextflow_args, load_project_configs
from cli.core.db import (
    DockedPose,
    DockingRun,
    Pose,
    ReceptorFrame,
    Substance,
    SubstanceState,
    Target,
    init_db,
)
from cli.core.runner import execute_nextflow

console = Console()


def _format_metric(val: float | str | None) -> str:
    if val is None or val == "" or str(val).lower() == "nan":
        return "-"
    f_val = float(val)
    if not math.isfinite(f_val):
        return "-"
    return f"{f_val:.2f}"


def _color_rmsd(val: float | None) -> str:
    if val is None or not math.isfinite(val):
        return "-"
    return (
        f"[bold green]{val:.2f}[/bold green]"
        if val <= 2.0
        else (
            f"[bold yellow]{val:.2f}[/bold yellow]"
            if val <= 3.0
            else f"[bold red]{val:.2f}[/bold red]"
        )
    )


def _resolve_frame_structural_poses(
    session, frame: ReceptorFrame
) -> tuple[Pose | None, list[Pose]]:
    """
    Deterministically resolves the reference pocket ligand pose and structural cofactor/ion poses
    strictly from structural origins ('crystal' or 'md_centroid'). Never mistakes docked poses for native.
    """
    ref_lig_pose = (
        session.query(Pose)
        .join(SubstanceState, Pose.substance_state_id == SubstanceState.id)
        .join(Substance, SubstanceState.substance_id == Substance.id)
        .filter(
            Pose.receptor_frame_id == frame.id,
            Pose.origin.in_(["crystal", "md_centroid"]),
            Substance.substance_type == "ligand",
        )
        .first()
    )
    if not ref_lig_pose and frame.parent_frame_id:
        ref_lig_pose = (
            session.query(Pose)
            .join(SubstanceState, Pose.substance_state_id == SubstanceState.id)
            .join(Substance, SubstanceState.substance_id == Substance.id)
            .filter(
                Pose.receptor_frame_id == frame.parent_frame_id,
                Pose.origin.in_(["crystal", "md_centroid"]),
                Substance.substance_type == "ligand",
            )
            .first()
        )

    cofactor_poses = (
        session.query(Pose)
        .join(SubstanceState, Pose.substance_state_id == SubstanceState.id)
        .join(Substance, SubstanceState.substance_id == Substance.id)
        .filter(
            Pose.receptor_frame_id == frame.id,
            Pose.origin.in_(["crystal", "md_centroid"]),
            Substance.substance_type.in_(["cofactor", "ion"]),
        )
        .all()
    )
    if not cofactor_poses and frame.parent_frame_id:
        cofactor_poses = (
            session.query(Pose)
            .join(SubstanceState, Pose.substance_state_id == SubstanceState.id)
            .join(Substance, SubstanceState.substance_id == Substance.id)
            .filter(
                Pose.receptor_frame_id == frame.parent_frame_id,
                Pose.origin.in_(["crystal", "md_centroid"]),
                Substance.substance_type.in_(["cofactor", "ion"]),
            )
            .all()
        )

    return ref_lig_pose, cofactor_poses


def _parse_mol_float_prop(mol: Chem.Mol, candidates: list[str]) -> float | None:
    prop_map = {p.strip().lower(): p for p in mol.GetPropNames()}
    for cand in candidates:
        c_low = cand.lower()
        if c_low in prop_map:
            raw = mol.GetProp(prop_map[c_low]).strip()
            if raw and raw.lower() != "nan":
                return float(raw)
    return None


def _ensure_substance_state_for_pose(
    session, cwd: Path, substance: Substance, mol: Chem.Mol
) -> SubstanceState:
    _, _, _, state_hash, isomeric_smiles, formal_charge = get_microstate_identity(mol)
    state_entity = (
        session.query(SubstanceState)
        .filter_by(substance_id=substance.id, state_hash=state_hash)
        .first()
    )
    if not state_entity:
        clean_res = substance.res_name or substance.name[:4].strip().upper()
        state_dir = cwd / "substances" / substance.inchikey / "states" / state_hash
        state_dir.mkdir(parents=True, exist_ok=True)
        state_sdf_file = state_dir / f"{clean_res}_state.sdf"
        with Chem.SDWriter(str(state_sdf_file)) as writer:
            writer.write(mol)

        state_entity = SubstanceState(
            substance_id=substance.id,
            state_hash=state_hash,
            isomeric_smiles=isomeric_smiles,
            formal_charge=formal_charge,
            state="raw",
            path_sdf=str(state_sdf_file.resolve()),
        )
        session.add(state_entity)
        session.flush()
    return state_entity


def _query_ready_structural_frames(
    session, target: str | None, frame_names: list[str] | None
) -> list[ReceptorFrame]:
    query = session.query(ReceptorFrame).filter(
        ReceptorFrame.state == "ready",
        ReceptorFrame.source_type.in_(
            ["crystal", "crystal_healed", "nmr", "md_centroid"]
        ),
    )
    if target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            session.close()
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif frame_names:
        query = query.filter(ReceptorFrame.name.in_(frame_names))

    frames = query.order_by(ReceptorFrame.name).all()
    if not frames:
        session.close()
        console.print(
            "[bold red]FATAL: No ready ReceptorFrames found. Run 'sbdd prep' first![/bold red]"
        )
        sys.exit(1)
    return frames


def _build_receptor_plan_entry(
    r: ReceptorFrame, ref_pose: Pose, cof_poses: list[Pose]
) -> dict:
    if not r.path_mae or not r.path_pdb:
        raise ValueError(
            f"FATAL: ReceptorFrame '{r.name}' is missing path_mae or path_pdb in 3NF database!"
        )
    path_mae = Path(r.path_mae).resolve()
    path_pdb = Path(r.path_pdb).resolve()
    if not path_mae.exists() or not path_pdb.exists():
        raise FileNotFoundError(
            f"FATAL: Receptor deliverables missing on disk for '{r.name}': {path_mae}, {path_pdb}"
        )

    if not ref_pose.path_sdf or not ref_pose.path_mae:
        raise ValueError(
            f"FATAL: Reference Pose '{ref_pose.name}' for '{r.name}' is missing path_sdf or path_mae!"
        )
    pose_sdf = Path(ref_pose.path_sdf).resolve()
    pose_mae = Path(ref_pose.path_mae).resolve()
    if not pose_sdf.exists() or not pose_mae.exists():
        raise FileNotFoundError(
            f"FATAL: Reference pose deliverables missing on disk: {pose_sdf}, {pose_mae}"
        )

    ion_sdfs, ion_maes, ion_resnames = [], [], []
    cof_sdfs, cof_maes, cof_resnames = [], [], []
    for c in cof_poses:
        if c.substance.substance_type == "ion":
            ion_sdfs.append(str(Path(c.path_sdf).resolve()))
            if c.path_mae:
                ion_maes.append(str(Path(c.path_mae).resolve()))
            ion_resnames.append(c.res_name)
        else:
            cof_sdfs.append(str(Path(c.path_sdf).resolve()))
            if c.path_mae:
                cof_maes.append(str(Path(c.path_mae).resolve()))
            cof_resnames.append(c.res_name)

    return {
        "name": r.name,
        "path_mae": str(path_mae),
        "path_pdb": str(path_pdb),
        "ref_sdf": str(pose_sdf),
        "ref_mae": str(pose_mae),
        "native_ligand": ref_pose.substance.name,
        "ion_sdfs": ion_sdfs,
        "ion_maes": ion_maes,
        "ion_resnames": ion_resnames,
        "cofactor_sdfs": cof_sdfs,
        "cofactor_maes": cof_maes,
        "cofactor_resnames": cof_resnames,
    }


def _ingest_docking_receipt(
    session,
    cwd: Path,
    run_name: str,
    mode: str,
    plan: dict,
    ref_mol_by_frame: dict[str, tuple[Substance, Chem.Mol]] | None = None,
) -> tuple[
    DockingRun,
    dict[tuple[str, str], dict[str, list[tuple[DockedPose, Chem.Mol]]]],
    int,
]:
    receipt_file = cwd / "dock_receipt.tsv"
    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not generated by Nextflow."
        )

    substance_by_title = {
        (
            item["ligand_title"] if "ligand_title" in item else item["name"]
        ): session.query(Substance).filter(Substance.id == item["substance_id"]).one()
        for item in plan["ligands"]
    }
    receipt_map: dict[tuple[str, str], Path] = {}
    with open(receipt_file, "r", encoding="utf-8") as rf:
        for line_num, line in enumerate(rf, start=1):
            line_clean = line.strip()
            if not line_clean:
                continue
            parts = line_clean.split("\t")
            if len(parts) < 3:
                session.close()
                raise ValueError(
                    f"FATAL: Corrupted line {line_num} in {receipt_file}: '{line_clean}'"
                )
            r_frame_id, r_engine, r_sdf_path = parts[0], parts[1], parts[2]
            sdf_file = Path(r_sdf_path).resolve()
            if not sdf_file.exists() or sdf_file.stat().st_size == 0:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: Published pose file specified in receipt does not exist or is empty: {sdf_file}"
                )
            receipt_map[(r_frame_id, r_engine)] = sdf_file

    dock_run = session.query(DockingRun).filter_by(name=run_name).first()
    if dock_run:
        session.delete(dock_run)
        session.flush()

    dock_run = DockingRun(name=run_name, mode=mode)
    session.add(dock_run)
    session.flush()

    parsed_modes: dict[
        tuple[str, str], dict[str, list[tuple[DockedPose, Chem.Mol]]]
    ] = {}
    total_modes_saved = 0

    for (frame_name, eng), sdf_file in sorted(receipt_map.items()):
        r_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
        if not r_frame:
            session.close()
            raise ValueError(
                f"FATAL: ReceptorFrame '{frame_name}' in receipt not found in database!"
            )

        frame_dir = sdf_file.parent
        glide_csv_scores: dict[str, list[float]] = {}
        glide_csv = frame_dir / f"{frame_name}_dock.csv"
        if eng == "glide" and glide_csv.exists() and glide_csv.stat().st_size > 0:
            with open(glide_csv, "r", encoding="utf-8", errors="ignore") as cf:
                for row in csv.DictReader(cf):
                    t_raw = (
                        row.get("title")
                        if row.get("title") is not None
                        else (
                            row.get("Title")
                            if row.get("Title") is not None
                            else row.get("s_m_title")
                        )
                    )
                    t_str = t_raw.strip() if t_raw else ""
                    sc_raw = (
                        row.get("r_i_docking_score")
                        if row.get("r_i_docking_score") is not None
                        else (
                            row.get("docking_score")
                            if row.get("docking_score") is not None
                            else row.get("r_i_glide_gscore")
                        )
                    )
                    if t_str and sc_raw:
                        if t_str not in glide_csv_scores:
                            glide_csv_scores[t_str] = []
                        glide_csv_scores[t_str].append(float(sc_raw))

        suppl = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False)
        mode_counters: dict[str, int] = {}

        for record_idx, m in enumerate(suppl, start=1):
            if m is None:
                session.close()
                raise ValueError(
                    f"FATAL: Corrupted molecule record at index {record_idx} in {sdf_file}!"
                )
            m_title = m.GetProp("_Name").strip() if m.HasProp("_Name") else ""
            if not m_title and m.HasProp("s_m_title"):
                m_title = m.GetProp("s_m_title").strip()
            if not m_title:
                session.close()
                raise ValueError(
                    f"FATAL: Unnamed molecule record at index {record_idx} in {sdf_file}!"
                )

            sub_obj = substance_by_title.get(m_title)
            if not sub_obj:
                session.close()
                raise KeyError(
                    f"FATAL: Molecule title '{m_title}' in {sdf_file} was not declared in docking plan!"
                )

            prev_count = mode_counters[m_title] if m_title in mode_counters else 0
            m_idx = prev_count + 1
            mode_counters[m_title] = m_idx

            g_score = _parse_mol_float_prop(
                m, ["r_i_docking_score", "docking_score", "r_i_glide_gscore", "gscore"]
            )
            if g_score is None and eng == "glide" and m_title in glide_csv_scores:
                csv_list = glide_csv_scores[m_title]
                if m_idx <= len(csv_list):
                    g_score = csv_list[m_idx - 1]

            v_aff = _parse_mol_float_prop(
                m, ["unidock_vina_score", "minimizedAffinity", "affinity"]
            )
            c_score = _parse_mol_float_prop(m, ["CNNscore"])
            c_aff = _parse_mol_float_prop(m, ["CNNaffinity"])

            mol_clean = Chem.Mol(m)
            mol_clean.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol_clean)

            insitu_rmsd_val: float | None = None
            if ref_mol_by_frame is not None and frame_name in ref_mol_by_frame:
                ref_pair = ref_mol_by_frame[frame_name]
                if ref_pair[0].id == sub_obj.id:
                    insitu_rmsd_val = round(
                        compute_symmetry_rmsd(ref_pair[1], mol_clean), 3
                    )

            dp = DockedPose(
                docking_run_id=dock_run.id,
                receptor_frame_id=r_frame.id,
                substance_id=sub_obj.id,
                substance_state_id=None,
                engine=eng,
                mode_index=m_idx,
                glide_score=g_score,
                vina_affinity=v_aff,
                cnn_score=c_score,
                cnn_affinity=c_aff,
                insitu_rmsd=insitu_rmsd_val,
                path_sdf=str(sdf_file),
            )
            session.add(dp)
            pair_key = (frame_name, m_title)
            if pair_key not in parsed_modes:
                parsed_modes[pair_key] = {}
            if eng not in parsed_modes[pair_key]:
                parsed_modes[pair_key][eng] = []
            parsed_modes[pair_key][eng].append((dp, mol_clean))
            total_modes_saved += 1

    for eng_dict in parsed_modes.values():
        glide_list = eng_dict["glide"] if "glide" in eng_dict else []
        unidock_list = eng_dict["unidock"] if "unidock" in eng_dict else []
        if glide_list and unidock_list:
            g_mols = [item[1] for item in glide_list]
            u_mols = [item[1] for item in unidock_list]
            min_for_g, min_for_u = compute_pairwise_symmetry_rmsd_matrix(g_mols, u_mols)
            for g_idx, (g_dp, _) in enumerate(glide_list):
                g_dp.cross_engine_rmsd = min_for_g[g_idx]
            for u_idx, (u_dp, _) in enumerate(unidock_list):
                u_dp.cross_engine_rmsd = min_for_u[u_idx]

    session.commit()
    return dock_run, parsed_modes, total_modes_saved


def _parse_csv_float(row: dict[str, str], key: str) -> float | None:
    if key not in row or row[key] is None:
        return None
    raw = str(row[key]).strip()
    if not raw or raw.lower() == "nan":
        return None
    return float(raw)


def _ingest_matrix_screening_results(
    session: Session,
    cwd: Path,
    run_name: str,
    receptors_plan: list[dict[str, str | list[str]]],
    ligands_plan: list[dict[str, str | int | bool]],
) -> tuple[DockingRun, list[dict[str, float | int | str | None]], int]:
    """
    Fast O(N) ingestion of virtual screening results directly from the worker-computed
    reports/ranking.csv and dock_receipt.tsv without re-parsing multi-gigabyte SDFs on master.
    """
    receipt_file = cwd / "dock_receipt.tsv"
    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not generated by Nextflow."
        )

    raw_csv = cwd / "reports" / "docking_combined_raw.csv"
    ranking_csv = cwd / "reports" / "ranking.csv"
    source_csv = raw_csv if (raw_csv.exists() and raw_csv.stat().st_size > 0) else ranking_csv
    if not source_csv.exists() or source_csv.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Master ranking table {source_csv} was not generated by SUMMARIZE_DOCKING worker."
        )

    receipt_map: dict[tuple[str, str], str] = {}
    with open(receipt_file, "r", encoding="utf-8") as rf:
        for line_num, line in enumerate(rf, start=1):
            line_clean = line.strip()
            if not line_clean:
                continue
            parts = line_clean.split("\t")
            if len(parts) < 3:
                session.close()
                raise ValueError(
                    f"FATAL: Corrupted line {line_num} in {receipt_file}: '{line_clean}'"
                )
            r_frame_id, r_engine, r_sdf_path = parts[0], parts[1], parts[2]
            sdf_file = Path(r_sdf_path).resolve()
            if not sdf_file.exists() or sdf_file.stat().st_size == 0:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: Published pose file specified in receipt does not exist or is empty: {sdf_file}"
                )
            receipt_map[(r_frame_id, r_engine)] = str(sdf_file)

    substance_id_by_title: dict[str, int] = {
        str(
            item["ligand_title"] if "ligand_title" in item else item["name"]
        ): int(item["substance_id"])
        for item in ligands_plan
    }
    frame_id_by_name: dict[str, int] = {}
    for rec in receptors_plan:
        f_name = str(rec["name"])
        r_frame = session.query(ReceptorFrame).filter_by(name=f_name).first()
        if not r_frame:
            session.close()
            raise ValueError(f"FATAL: ReceptorFrame '{f_name}' not found in database!")
        frame_id_by_name[f_name] = r_frame.id

    dock_run = session.query(DockingRun).filter_by(name=run_name).first()
    if dock_run:
        session.delete(dock_run)
        session.flush()

    dock_run = DockingRun(name=run_name, mode="matrix")
    session.add(dock_run)
    session.flush()

    summary_rows: list[dict[str, float | int | str | None]] = []
    dp_objects: list[DockedPose] = []

    with open(source_csv, "r", encoding="utf-8") as cf:
        for row in csv.DictReader(cf):
            m_title = str(row["Title"]).strip() if "Title" in row and row["Title"] else ""
            frame_name = (
                str(row["Receptor"]).strip()
                if "Receptor" in row and row["Receptor"]
                else (str(row["Frame"]).strip() if "Frame" in row and row["Frame"] else "")
            )
            if not m_title or not frame_name:
                continue
            if m_title not in substance_id_by_title:
                session.close()
                raise KeyError(
                    f"FATAL: Substance '{m_title}' in {source_csv} was not declared in docking plan!"
                )
            if frame_name not in frame_id_by_name:
                session.close()
                raise KeyError(
                    f"FATAL: ReceptorFrame '{frame_name}' in {source_csv} was not declared in docking plan!"
                )

            sub_id = substance_id_by_title[m_title]
            f_id = frame_id_by_name[frame_name]

            g_score = _parse_csv_float(row, "GlideScore")
            g_cnn_sc = _parse_csv_float(row, "Glide_CNNscore")
            g_cnn_aff = _parse_csv_float(row, "Glide_CNNaffinity")
            g_modes_f = _parse_csv_float(row, "Glide_Modes")
            g_best_m = _parse_csv_float(row, "Glide_Best_Mode")
            n_glide = int(g_modes_f) if g_modes_f is not None else (1 if g_score is not None else 0)
            g_mode_idx = int(g_best_m) if g_best_m is not None else 1

            u_vina = _parse_csv_float(row, "UniDock_Vina")
            u_cnn_sc = _parse_csv_float(row, "UniDock_CNNscore")
            u_cnn_aff = _parse_csv_float(row, "UniDock_CNNaffinity")
            u_modes_f = _parse_csv_float(row, "UniDock_Modes")
            u_best_m = _parse_csv_float(row, "UniDock_Best_Mode")
            n_unidock = int(u_modes_f) if u_modes_f is not None else (1 if u_vina is not None else 0)
            u_mode_idx = int(u_best_m) if u_best_m is not None else 1

            c_size_f = _parse_csv_float(row, "Cluster_Size")
            n_clust_f = _parse_csv_float(row, "Num_Clusters")
            c_size = int(c_size_f) if c_size_f is not None else 1
            n_clust = int(n_clust_f) if n_clust_f is not None else 1

            pair_rmsd = _parse_csv_float(row, "Pose_RMSD")
            ult_score = _parse_csv_float(row, "Ultimate_Score")

            if g_score is not None and (frame_name, "glide") in receipt_map:
                dp_objects.append(
                    DockedPose(
                        docking_run_id=dock_run.id,
                        receptor_frame_id=f_id,
                        substance_id=sub_id,
                        substance_state_id=None,
                        engine="glide",
                        mode_index=g_mode_idx,
                        glide_score=g_score,
                        vina_affinity=None,
                        cnn_score=g_cnn_sc,
                        cnn_affinity=g_cnn_aff,
                        cross_engine_rmsd=pair_rmsd,
                        consensus_score=ult_score,
                        path_sdf=receipt_map[(frame_name, "glide")],
                    )
                )

            if u_vina is not None and (frame_name, "unidock") in receipt_map:
                dp_objects.append(
                    DockedPose(
                        docking_run_id=dock_run.id,
                        receptor_frame_id=f_id,
                        substance_id=sub_id,
                        substance_state_id=None,
                        engine="unidock",
                        mode_index=u_mode_idx,
                        glide_score=None,
                        vina_affinity=u_vina,
                        cnn_score=u_cnn_sc,
                        cnn_affinity=u_cnn_aff,
                        cross_engine_rmsd=pair_rmsd,
                        consensus_score=ult_score,
                        path_sdf=receipt_map[(frame_name, "unidock")],
                    )
                )

            summary_rows.append(
                {
                    "frame": frame_name,
                    "title": m_title,
                    "n_glide": n_glide,
                    "n_unidock": n_unidock,
                    "g_mode": g_mode_idx,
                    "u_mode": u_mode_idx,
                    "c_size": c_size,
                    "n_clust": n_clust,
                    "g_score": g_score,
                    "g_cnn_aff": g_cnn_aff,
                    "g_cnn_sc": g_cnn_sc,
                    "u_vina": u_vina,
                    "u_cnn_aff": u_cnn_aff,
                    "u_cnn_sc": u_cnn_sc,
                    "pair_rmsd": pair_rmsd,
                    "ult_score": ult_score,
                }
            )

    session.bulk_save_objects(dp_objects)
    session.commit()
    return dock_run, summary_rows, len(dp_objects)


def run_dock(
    target: str | None = None,
    frame_names: list[str] | None = None,
    batch_names: list[str] | None = None,
    ligand_names: list[str] | None = None,
    name: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """
    Runs Virtual Screening / Cross-Docking of library substances into ReceptorFrames.
    Selects representative poses strictly by best scoring energy/affinity and stores all
    modes (1..16) in 3NF DockingRun / DockedPose tables without native RMSD columns.
    """
    cwd = Path.cwd()
    session = init_db()

    frames = _query_ready_structural_frames(session, target, frame_names)

    lig_query = session.query(Substance).filter_by(substance_type="ligand", state="ready")
    if batch_names:
        lig_query = lig_query.filter(Substance.batch_id.in_(batch_names))

    if not ligand_names or ligand_names == ["all"]:
        ligands = lig_query.order_by(Substance.name).all()
        if not ligands:
            scope_desc = f" for batch(es) '{batch_names}'" if batch_names else ""
            session.close()
            console.print(
                f"[bold red]FATAL: No ready ligands found in database{scope_desc}! Run 'sbdd add ligand' and 'sbdd prep' first.[/bold red]"
            )
            sys.exit(1)
    else:
        ligands = lig_query.filter(Substance.name.in_(ligand_names)).order_by(Substance.name).all()
        if len(ligands) != len(set(ligand_names)):
            found = {l.name for l in ligands}
            missing = set(ligand_names) - found
            session.close()
            console.print(
                f"[bold red]FATAL: Ligands not found in DB or not ready: {missing}[/bold red]"
            )
            sys.exit(1)

    plan = {
        "mode": "matrix",
        "receptors": [],
        "ligands": [],
        "pairs": [],
    }

    for r in frames:
        pose_obj, cof_poses = _resolve_frame_structural_poses(session, r)
        if not pose_obj or not pose_obj.substance:
            session.close()
            console.print(
                f"[bold red]FATAL: ReceptorFrame '{r.name}' has no structural reference pose ('crystal' or 'md_centroid') to define the binding pocket![/bold red]"
            )
            sys.exit(1)
        plan["receptors"].append(_build_receptor_plan_entry(r, pose_obj, cof_poses))

    for l in ligands:
        if l.state != "ready" or not l.path_mae or not l.path_sdf:
            session.close()
            console.print(
                f"[bold red]FATAL: Chemical Ligand '{l.name}' is in state '{l.state}'! Run 'sbdd prep' first.[/bold red]"
            )
            sys.exit(1)

        lig_mae = Path(l.path_mae).resolve()
        lig_sdf = Path(l.path_sdf).resolve()
        if not lig_mae.exists() or not lig_sdf.exists():
            session.close()
            raise FileNotFoundError(
                f"FATAL: Prepared files missing for '{l.name}' ({lig_mae}). Run 'sbdd prep' first."
            )

        plan["ligands"].append(
            {
                "name": l.name,
                "path_sdf": str(lig_sdf),
                "path_mae": str(lig_mae),
                "is_native": False,
                "ligand_title": l.name,
                "substance_id": l.id,
            }
        )

    run_name = name.strip() if name else f"screen_{len(frames)}x{len(ligands)}"
    console.print(
        f"[bold cyan]Cross-Docking Matrix Plan generated: {len(frames)} Receptor Frame(s) × {len(ligands)} Stereospecific Substance(s)[/bold cyan]"
    )

    plan_file = cwd / "dock_plan.json"
    with open(plan_file, "w", encoding="utf-8") as f:
        json.dump(plan, f, indent=2)
    (cwd / "dock_receipt.tsv").unlink(missing_ok=True)

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="screen")
    nf_args.extend(
        [
            "--dock_plan",
            str(plan_file.resolve()),
            "--dock_run_name",
            run_name,
        ]
    )
    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
    ret_code = execute_nextflow(
        "screen", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )
    if ret_code != 0:
        session.close()
        console.print("[bold red]Nextflow docking failed.[/bold red]")
        sys.exit(ret_code)

    _, summary_rows, total_modes_saved = _ingest_matrix_screening_results(
        session=session,
        cwd=cwd,
        run_name=run_name,
        receptors_plan=plan["receptors"],
        ligands_plan=plan["ligands"],
    )

    has_glide = any(r["g_score"] is not None for r in summary_rows)
    has_unidock = any(r["u_vina"] is not None for r in summary_rows)
    dual_engine = has_glide and has_unidock

    # 1. Render Campaign Overview per Receptor Frame
    overview_table = Table(
        title=f"📊 Virtual Screening Campaign Overview ('{run_name}')",
        show_header=True,
    )
    overview_table.add_column("Frame", style="bold cyan")
    if dual_engine:
        overview_table.add_column("Docked (Glide / UniDock)", justify="center")
        overview_table.add_column("Consensus Poses (ΔRMSD ≤ 2.0 Å)", justify="center", style="bold green")
        overview_table.add_column("GlideScore (Best / Median)", justify="right")
        overview_table.add_column("UniDock Vina (Best / Median)", justify="right")
        overview_table.add_column("CNNscore (Best / Median)", justify="right", style="magenta")
    elif has_unidock:
        overview_table.add_column("Substances Docked", justify="center")
        overview_table.add_column("Converged Basins (≥ 4 poses)", justify="center", style="bold green")
        overview_table.add_column("UniDock Vina (Best / Median)", justify="right")
        overview_table.add_column("CNNscore (Best / Median)", justify="right", style="magenta")
        overview_table.add_column("CNNaffinity (Best / Median)", justify="right")
    else:
        overview_table.add_column("Substances Docked", justify="center")
        overview_table.add_column("Converged Basins (≥ 4 poses)", justify="center", style="bold green")
        overview_table.add_column("GlideScore (Best / Median)", justify="right")
        overview_table.add_column("CNNscore (Best / Median)", justify="right", style="magenta")

    for rec in plan["receptors"]:
        f_name = rec["name"]
        f_rows = [r for r in summary_rows if r["frame"] == f_name]
        g_vals = [float(r["g_score"]) for r in f_rows if r["g_score"] is not None]
        u_vals = [float(r["u_vina"]) for r in f_rows if r["u_vina"] is not None]
        rmsd_vals = [float(r["pair_rmsd"]) for r in f_rows if r["pair_rmsd"] is not None]
        cnn_sc_vals = [
            float(v)
            for r in f_rows
            for v in (r["g_cnn_sc"], r["u_cnn_sc"])
            if v is not None
        ]
        cnn_aff_vals = [
            float(v)
            for r in f_rows
            for v in (r["g_cnn_aff"], r["u_cnn_aff"])
            if v is not None
        ]

        g_stat_str = (
            f"{min(g_vals):.2f} / {statistics.median(g_vals):.2f}" if g_vals else "-"
        )
        u_stat_str = (
            f"{min(u_vals):.2f} / {statistics.median(u_vals):.2f}" if u_vals else "-"
        )
        cnn_stat_str = (
            f"{max(cnn_sc_vals):.2f} / {statistics.median(cnn_sc_vals):.2f}"
            if cnn_sc_vals
            else "-"
        )
        aff_stat_str = (
            f"{max(cnn_aff_vals):.2f} / {statistics.median(cnn_aff_vals):.2f}"
            if cnn_aff_vals
            else "-"
        )

        if dual_engine:
            consensus_cnt = sum(1 for v in rmsd_vals if v <= 2.0)
            consensus_pct = (100.0 * consensus_cnt / len(rmsd_vals)) if rmsd_vals else 0.0
            consensus_str = (
                f"{consensus_cnt}/{len(rmsd_vals)} ({consensus_pct:.1f}%)"
                if rmsd_vals
                else "-"
            )
            overview_table.add_row(
                f_name,
                f"{len(g_vals)} / {len(u_vals)}",
                consensus_str,
                g_stat_str,
                u_stat_str,
                cnn_stat_str,
            )
        elif has_unidock:
            dense_cnt = sum(1 for r in f_rows if int(r["c_size"]) >= 4)
            dense_pct = (100.0 * dense_cnt / len(f_rows)) if f_rows else 0.0
            overview_table.add_row(
                f_name,
                str(len(u_vals)),
                f"{dense_cnt}/{len(f_rows)} ({dense_pct:.1f}%)",
                u_stat_str,
                cnn_stat_str,
                aff_stat_str,
            )
        else:
            dense_cnt = sum(1 for r in f_rows if int(r["c_size"]) >= 4)
            dense_pct = (100.0 * dense_cnt / len(f_rows)) if f_rows else 0.0
            overview_table.add_row(
                f_name,
                str(len(g_vals)),
                f"{dense_cnt}/{len(f_rows)} ({dense_pct:.1f}%)",
                g_stat_str,
                cnn_stat_str,
            )

    # 2. Deduplicate across frames: select the single best conformational frame per unique substance
    best_by_title: dict[str, dict[str, float | int | str | None]] = {}
    for r in summary_rows:
        t_str = str(r["title"])
        u_val = float(r["ult_score"]) if r["ult_score"] is not None else -1.0
        if t_str not in best_by_title:
            best_by_title[t_str] = r
        else:
            prev_u = (
                float(best_by_title[t_str]["ult_score"])
                if best_by_title[t_str]["ult_score"] is not None
                else -1.0
            )
            if u_val > prev_u:
                best_by_title[t_str] = r

    dedup_rows = sorted(
        list(best_by_title.values()),
        key=lambda x: float(x["ult_score"]) if x["ult_score"] is not None else -1.0,
        reverse=True,
    )

    top_k = min(25, len(dedup_rows))
    table = Table(
        title=f"🎯 Top-{top_k} Unique Substances in Best Frame (Ensemble Clustered | Total DB Records: {total_modes_saved})",
        show_header=True,
    )
    table.add_column("Rank", justify="right", style="bold")
    table.add_column("Best Frame", style="cyan")
    table.add_column("Substance", style="bold")

    if dual_engine:
        table.add_column("Basin Modes (G/U)", justify="center", style="dim")
        table.add_column("Glide Score", justify="right")
        table.add_column("Glide CNNscore (aff)", justify="right")
        table.add_column("UniDock Vina", justify="right")
        table.add_column("UniDock CNNscore (aff)", justify="right")
        table.add_column("Basin ΔRMSD (Å)", justify="center")
        table.add_column("Ultimate Score", justify="right", style="bold magenta")

        for rank_idx, r in enumerate(dedup_rows[:top_k], start=1):
            g_cnn_str = (
                f"{float(r['g_cnn_sc']):.2f} ({float(r['g_cnn_aff']):.2f})"
                if (r["g_cnn_sc"] is not None and r["g_cnn_aff"] is not None)
                else _format_metric(r["g_cnn_sc"])
            )
            u_cnn_str = (
                f"{float(r['u_cnn_sc']):.2f} ({float(r['u_cnn_aff']):.2f})"
                if (r["u_cnn_sc"] is not None and r["u_cnn_aff"] is not None)
                else _format_metric(r["u_cnn_sc"])
            )
            pair_rmsd_val = float(r["pair_rmsd"]) if r["pair_rmsd"] is not None else None
            modes_str = f"#{r['g_mode']}/#{r['u_mode']} ({r['c_size']}p/{r['n_clust']}c)"
            table.add_row(
                str(rank_idx),
                str(r["frame"]),
                str(r["title"]),
                modes_str,
                _format_metric(r["g_score"]),
                g_cnn_str,
                _format_metric(r["u_vina"]),
                u_cnn_str,
                _color_rmsd(pair_rmsd_val),
                _format_metric(r["ult_score"]),
            )
    elif has_unidock:
        table.add_column("Basin Mode", justify="center", style="dim")
        table.add_column("UniDock Vina (kcal/mol)", justify="right")
        table.add_column("CNNscore (Validity)", justify="right", style="green")
        table.add_column("CNNaffinity (pKd)", justify="right")
        table.add_column("Ultimate Score", justify="right", style="bold magenta")

        for rank_idx, r in enumerate(dedup_rows[:top_k], start=1):
            modes_str = f"#{r['u_mode']} ({r['c_size']}p/{r['n_clust']}c)"
            table.add_row(
                str(rank_idx),
                str(r["frame"]),
                str(r["title"]),
                modes_str,
                _format_metric(r["u_vina"]),
                _format_metric(r["u_cnn_sc"]),
                _format_metric(r["u_cnn_aff"]),
                _format_metric(r["ult_score"]),
            )
    else:
        table.add_column("Basin Mode", justify="center", style="dim")
        table.add_column("Glide Score (kcal/mol)", justify="right")
        table.add_column("CNNscore (Validity)", justify="right", style="green")
        table.add_column("CNNaffinity (pKd)", justify="right")
        table.add_column("Ultimate Score", justify="right", style="bold magenta")

        for rank_idx, r in enumerate(dedup_rows[:top_k], start=1):
            modes_str = f"#{r['g_mode']} ({r['c_size']}p/{r['n_clust']}c)"
            table.add_row(
                str(rank_idx),
                str(r["frame"]),
                str(r["title"]),
                modes_str,
                _format_metric(r["g_score"]),
                _format_metric(r["g_cnn_sc"]),
                _format_metric(r["g_cnn_aff"]),
                _format_metric(r["ult_score"]),
            )

    session.close()
    console.print(
        f"\n[bold green]✓ Virtual Screening Run '{run_name}' saved to 3NF DB ({total_modes_saved} records across {len(summary_rows)} pair(s), 0 Systems spammed).[/bold green]\n"
    )
    console.print(overview_table)
    console.print()
    console.print(table)
    console.print(
        f"\n[dim]📄 Full ranked table ({len(summary_rows)} rows): {cwd / 'reports' / 'ranking.csv'}[/dim]"
    )


def run_redock(
    target: str | None = None,
    frame_names: list[str] | None = None,
    name: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """
    Runs 1-to-1 Redocking of each ReceptorFrame's native structural ligand ('crystal' or 'md_centroid')
    back into its own pocket. Evaluates symmetry-aware in-situ RMSD via compute_symmetry_rmsd
    across all 16 generated modes and selects the best-reproducing pose per engine.
    """
    cwd = Path.cwd()
    session = init_db()

    frames = _query_ready_structural_frames(session, target, frame_names)

    active_frames = []
    for r in frames:
        native_pose, cof_poses = _resolve_frame_structural_poses(session, r)
        if not native_pose or not native_pose.substance:
            continue
        active_frames.append((r, native_pose, native_pose.substance, cof_poses))

    if not active_frames:
        session.close()
        console.print(
            "[bold red]FATAL: None of the selected ready frames have a native structural ligand pose ('crystal' or 'md_centroid') for redocking![/bold red]"
        )
        sys.exit(1)

    plan = {
        "mode": "self",
        "receptors": [],
        "ligands": [],
        "pairs": [],
    }
    ref_mol_by_frame: dict[str, tuple[Substance, Chem.Mol]] = {}

    for r, native_pose, native_lig, cof_poses in active_frames:
        plan["receptors"].append(_build_receptor_plan_entry(r, native_pose, cof_poses))

        pose_sdf = Path(native_pose.path_sdf).resolve()
        ref_suppl = Chem.SDMolSupplier(str(pose_sdf), removeHs=False, sanitize=False)
        ref_mol = next((m for m in ref_suppl if m is not None), None)
        if ref_mol is None:
            session.close()
            raise ValueError(f"FATAL: Failed to parse reference SDF {pose_sdf}")
        ref_mol_by_frame[r.name] = (native_lig, ref_mol)

        if (
            native_lig.state != "ready"
            or not native_lig.path_mae
            or not native_lig.path_sdf
        ):
            session.close()
            console.print(
                f"[bold red]FATAL: Chemical Substance '{native_lig.name}' is in state '{native_lig.state}'! "
                "Run 'sbdd prep' first so LigPrep generates canonical 3D states.[/bold red]"
            )
            sys.exit(1)

        lig_mae = Path(native_lig.path_mae).resolve()
        lig_sdf = Path(native_lig.path_sdf).resolve()
        if not lig_mae.exists() or not lig_sdf.exists():
            session.close()
            raise FileNotFoundError(
                f"FATAL: LigPrep deliverables missing for '{native_lig.name}': {lig_mae}"
            )

        plan["ligands"].append(
            {
                "name": f"{r.name}_{native_lig.name}",
                "path_sdf": str(lig_sdf),
                "path_mae": str(lig_mae),
                "is_native": True,
                "receptor": r.name,
                "ligand_title": native_lig.name,
                "substance_id": native_lig.id,
            }
        )

    run_name = name.strip() if name else f"redock_{len(active_frames)}"
    console.print(
        f"[bold green]✓ 1-to-1 Redocking Mode: Matched {len(active_frames)} receptor frame(s) directly to their native ligands.[/bold green]"
    )

    plan_file = cwd / "dock_plan.json"
    with open(plan_file, "w", encoding="utf-8") as f:
        json.dump(plan, f, indent=2)
    (cwd / "dock_receipt.tsv").unlink(missing_ok=True)

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="screen")
    nf_args.extend(
        [
            "--dock_plan",
            str(plan_file.resolve()),
            "--dock_run_name",
            run_name,
        ]
    )
    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
    ret_code = execute_nextflow(
        "screen", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )
    if ret_code != 0:
        session.close()
        console.print("[bold red]Nextflow redocking failed.[/bold red]")
        sys.exit(ret_code)

    _, parsed_modes, total_modes_saved = _ingest_docking_receipt(
        session=session,
        cwd=cwd,
        run_name=run_name,
        mode="self",
        plan=plan,
        ref_mol_by_frame=ref_mol_by_frame,
    )

    summary_rows = []
    for (frame_name, m_title), eng_dict in parsed_modes.items():
        glide_list = eng_dict.get("glide", [])
        unidock_list = eng_dict.get("unidock", [])

        # Select representative pose strictly by minimum symmetry-aware In-Situ RMSD to native pose
        g_best = (
            min(
                glide_list,
                key=lambda x: (
                    x[0].insitu_rmsd if x[0].insitu_rmsd is not None else float("inf")
                ),
            )
            if glide_list
            else None
        )
        u_best = (
            min(
                unidock_list,
                key=lambda x: (
                    x[0].insitu_rmsd if x[0].insitu_rmsd is not None else float("inf")
                ),
            )
            if unidock_list
            else None
        )

        g_top = g_best[0] if g_best else None
        u_top = u_best[0] if u_best else None

        best_pair_rmsd = None
        if g_best and u_best:
            best_pair_rmsd = round(compute_symmetry_rmsd(g_best[1], u_best[1]), 3)

        summary_rows.append(
            {
                "frame": frame_name,
                "title": m_title,
                "n_glide": len(glide_list),
                "n_unidock": len(unidock_list),
                "g_mode": g_top.mode_index if g_top else None,
                "g_score": g_top.glide_score if g_top else None,
                "g_cnn_aff": g_top.cnn_affinity if g_top else None,
                "g_cnn_sc": g_top.cnn_score if g_top else None,
                "g_rmsd": g_top.insitu_rmsd if g_top else None,
                "u_mode": u_top.mode_index if u_top else None,
                "u_vina": u_top.vina_affinity if u_top else None,
                "u_cnn_aff": u_top.cnn_affinity if u_top else None,
                "u_cnn_sc": u_top.cnn_score if u_top else None,
                "u_rmsd": u_top.insitu_rmsd if u_top else None,
                "pair_rmsd": best_pair_rmsd,
            }
        )

    summary_rows.sort(
        key=lambda r: min(
            r["g_rmsd"] if r["g_rmsd"] is not None else float("inf"),
            r["u_rmsd"] if r["u_rmsd"] is not None else float("inf"),
        )
    )

    table = Table(
        title=f"🎯 Redocking Validation Run '{run_name}' Summary (Best In-Situ RMSD Poses | Total Modes in DB: {total_modes_saved})",
        show_header=True,
    )
    table.add_column("Frame", style="cyan")
    table.add_column("Native Ligand", style="bold")
    table.add_column("Glide Best Mode", justify="center")
    table.add_column("Glide Score", justify="right")
    table.add_column("Glide CNNaff (p)", justify="right")
    table.add_column("Glide Ref RMSD (Å)", justify="center")
    table.add_column("UniDock Best Mode", justify="center")
    table.add_column("UniDock Vina", justify="right")
    table.add_column("UniDock CNNaff (p)", justify="right")
    table.add_column("UniDock Ref RMSD (Å)", justify="center")
    table.add_column("G↔U ΔRMSD (Å)", justify="center")
    table.add_column("Verdict", justify="center")

    for r in summary_rows:
        g_mode_str = f"#{r['g_mode']}/{r['n_glide']}" if r["g_mode"] else "-"
        u_mode_str = f"#{r['u_mode']}/{r['n_unidock']}" if r["u_mode"] else "-"
        g_cnn_str = (
            f"{r['g_cnn_aff']:.2f} ({r['g_cnn_sc']:.2f})"
            if (r["g_cnn_aff"] is not None and r["g_cnn_sc"] is not None)
            else _format_metric(r["g_cnn_aff"])
        )
        u_cnn_str = (
            f"{r['u_cnn_aff']:.2f} ({r['u_cnn_sc']:.2f})"
            if (r["u_cnn_aff"] is not None and r["u_cnn_sc"] is not None)
            else _format_metric(r["u_cnn_aff"])
        )

        vg = r["g_rmsd"]
        vu = r["u_rmsd"]
        best_r = min(
            vg if vg is not None else float("inf"),
            vu if vu is not None else float("inf"),
        )
        if best_r <= 2.0:
            verdict = "[bold green]✓ SUCCESS[/bold green]"
        elif math.isfinite(best_r):
            verdict = "[bold yellow]⚠ DRIFT[/bold yellow]"
        else:
            verdict = "[dim]-[/dim]"

        table.add_row(
            r["frame"],
            r["title"],
            g_mode_str,
            _format_metric(r["g_score"]),
            g_cnn_str,
            _color_rmsd(vg),
            u_mode_str,
            _format_metric(r["u_vina"]),
            u_cnn_str,
            _color_rmsd(vu),
            _color_rmsd(r["pair_rmsd"]),
            verdict,
        )

    session.close()
    console.print(
        f"\n[bold green]✓ Redocking Run '{run_name}' saved to 3NF DB ({total_modes_saved} modes across {len(summary_rows)} frame(s), 0 Systems spammed).[/bold green]\n"
    )
    console.print(table)
