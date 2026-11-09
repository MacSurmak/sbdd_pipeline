"""
File: cli/commands/dock.py
Description: CLI handler for 'sbdd dock'. Reads DB to construct Cross-Docking plan.
"""

import csv
import json
import math
from pathlib import Path
import sys
from rdkit import Chem
from rich.console import Console
from rich.table import Table
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
                mol_block=Chem.MolToMolBlock(mol_clean),
            )
            session.add(dp)
            pair_key = (frame_name, m_title)
            if pair_key not in parsed_modes:
                parsed_modes[pair_key] = {}
            if eng not in parsed_modes[pair_key]:
                parsed_modes[pair_key][eng] = []
            parsed_modes[pair_key][eng].append((dp, mol_clean))
            total_modes_saved += 1

    # Compute exact symmetry-aware Cross-Engine RMSD between Glide modes and Uni-Dock modes
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

    _, parsed_modes, total_modes_saved = _ingest_docking_receipt(
        session=session,
        cwd=cwd,
        run_name=run_name,
        mode="matrix",
        plan=plan,
        ref_mol_by_frame=None,
    )

    summary_rows = []
    for (frame_name, m_title), eng_dict in parsed_modes.items():
        glide_list = eng_dict.get("glide", [])
        unidock_list = eng_dict.get("unidock", [])

        # Select representative pose strictly by each engine's native physical score
        g_best = (
            min(
                glide_list,
                key=lambda x: (
                    x[0].glide_score if x[0].glide_score is not None else float("inf")
                ),
            )
            if glide_list
            else None
        )
        u_best = (
            min(
                unidock_list,
                key=lambda x: (
                    x[0].vina_affinity
                    if x[0].vina_affinity is not None
                    else float("inf")
                ),
            )
            if unidock_list
            else None
        )

        g_top = g_best[0] if g_best else None
        u_top = u_best[0] if u_best else None

        top_pair_rmsd = None
        if g_best and u_best:
            top_pair_rmsd = round(compute_symmetry_rmsd(g_best[1], u_best[1]), 3)

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
                "u_mode": u_top.mode_index if u_top else None,
                "u_vina": u_top.vina_affinity if u_top else None,
                "u_cnn_aff": u_top.cnn_affinity if u_top else None,
                "u_cnn_sc": u_top.cnn_score if u_top else None,
                "pair_rmsd": top_pair_rmsd,
            }
        )

    summary_rows.sort(
        key=lambda r: (
            (r["g_score"] if r["g_score"] is not None else 0.0)
            + (r["u_vina"] if r["u_vina"] is not None else 0.0)
        )
    )

    table = Table(
        title=f"🎯 Virtual Screening Run '{run_name}' Summary (Best Native Score Poses | Total Modes in DB: {total_modes_saved})",
        show_header=True,
    )
    table.add_column("Rank", justify="right", style="bold")
    table.add_column("Frame", style="cyan")
    table.add_column("Substance", style="bold")
    table.add_column("Modes (G/U)", justify="center", style="dim")
    table.add_column("Glide Score", justify="right")
    table.add_column("Glide CNNaff (p)", justify="right")
    table.add_column("UniDock Vina", justify="right")
    table.add_column("UniDock CNNaff (p)", justify="right")
    table.add_column("Glide↔UniDock ΔRMSD (Å)", justify="center")

    for rank_idx, r in enumerate(summary_rows, start=1):
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
        table.add_row(
            str(rank_idx),
            r["frame"],
            r["title"],
            f"{r['n_glide']}/{r['n_unidock']}",
            _format_metric(r["g_score"]),
            g_cnn_str,
            _format_metric(r["u_vina"]),
            u_cnn_str,
            _color_rmsd(r["pair_rmsd"]),
        )

    session.close()
    console.print(
        f"\n[bold green]✓ Virtual Screening Run '{run_name}' saved to 3NF DB ({total_modes_saved} modes across {len(summary_rows)} pair(s), 0 Systems spammed).[/bold green]\n"
    )
    console.print(table)


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
    table.add_column("Glide↔UniDock ΔRMSD (Å)", justify="center")
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
