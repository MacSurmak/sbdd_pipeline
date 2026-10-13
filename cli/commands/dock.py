"""
File: cli/commands/dock.py
Description: CLI handler for 'sbdd dock'. Reads DB to construct Cross-Docking plan.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def _format_metric(val):
    if val is None or val == "" or str(val).lower() == "nan":
        return "-"
    try:
        return f"{float(val):.2f}"
    except (ValueError, TypeError):
        return str(val)


def _color_rmsd(val_str):
    if val_str is None or val_str == "" or str(val_str).lower() in ["nan", "-"]:
        return "-"
    try:
        v = float(val_str)
        return (
            f"[bold green]{v:.2f}[/bold green]"
            if v <= 2.0
            else f"[bold red]{v:.2f}[/bold red]"
        )
    except (ValueError, TypeError):
        return str(val_str)


def run_dock(
    target: str | None = None,
    frame_names: list[str] | None = None,
    ligand_names: list[str] | None = None,
    self_dock: bool = False,
    name: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """Run distributed docking into ReceptorFrames."""
    cwd = Path.cwd()
    session = init_db()
    is_self = self_dock

    # 1. Fetch requested Ready Frames
    query = session.query(ReceptorFrame).filter_by(state="ready")
    if target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif frame_names:
        query = query.filter(ReceptorFrame.name.in_(frame_names))

    frames = query.all()
    if not frames:
        console.print(
            "[bold red]FATAL: No ready ReceptorFrames found to dock. Run 'sbdd prep' first![/bold red]"
        )
        sys.exit(1)

    # 2. Build Plan (Self 1-to-1 vs Matrix Screening)
    plan = {
        "mode": "self" if is_self else "matrix",
        "receptors": [],
        "ligands": [],
        "pairs": [],
    }

    if is_self:
        active_frames = []

        for r in frames:
            target_sys = (
                session.query(System)
                .filter(System.receptor_frame_id == r.id, System.pose_id.isnot(None))
                .first()
            )
            if not target_sys or not target_sys.pose:
                continue

            native_pose = target_sys.pose
            native_lig = native_pose.ligand

            active_frames.append((r, native_pose, native_lig))

        if not active_frames:
            console.print(
                f"[bold green]✓ All {len(frames)} ready frame(s) already have redocking poses completed. Nothing to do![/bold green]"
            )
            session.close()
            return

        for r, native_pose, native_lig in active_frames:
            # Pure 3NF Asset Resolution: zero heuristics, zero with_suffix()
            if not r.path_mae or not r.path_pdb:
                raise ValueError(
                    f"FATAL: ReceptorFrame '{r.name}' does not have both path_mae and path_pdb "
                    f"registered in the 3NF database! Entity must be prepared cleanly."
                )

            path_mae = Path(r.path_mae).resolve()
            path_pdb = Path(r.path_pdb).resolve()

            if not path_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered receptor MAE does not exist: {path_mae}"
                )
            if not path_pdb.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered receptor PDB does not exist: {path_pdb}"
                )

            if not native_pose.path_sdf or not native_pose.path_mae:
                raise ValueError(
                    f"FATAL: Native Pose for '{r.name}' does not have both path_sdf and path_mae "
                    f"registered in the 3NF database!"
                )

            pose_sdf = Path(native_pose.path_sdf).resolve()
            pose_mae = Path(native_pose.path_mae).resolve()

            if not pose_sdf.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered pose SDF does not exist: {pose_sdf}"
                )
            if not pose_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered pose MAE does not exist: {pose_mae}"
                )

            plan["receptors"].append(
                {
                    "name": r.name,
                    "path_mae": str(path_mae),
                    "path_pdb": str(path_pdb),
                    "ref_sdf": str(pose_sdf),
                    "ref_mae": str(pose_mae),
                    "native_ligand": native_lig.name,
                }
            )

            unique_lig_name = f"{r.name}_{native_lig.name}"

            # Strict validation: Candidate ligand MUST be in state 'ready' (processed via LigPrep)
            if (
                native_lig.state != "ready"
                or not native_lig.path_mae
                or not native_lig.path_sdf
            ):
                console.print(
                    f"[bold red]FATAL: Chemical Ligand '{native_lig.name}' is in state '{native_lig.state}'! "
                    f"Run 'sbdd prep' first so LigPrep generates canonical 3D states.[/bold red]"
                )
                session.close()
                sys.exit(1)

            lig_mae = Path(native_lig.path_mae).resolve()
            lig_sdf = Path(native_lig.path_sdf).resolve()

            if not lig_mae.exists() or not lig_sdf.exists():
                console.print(
                    f"[bold red]FATAL: LigPrep deliverables missing for '{native_lig.name}': {lig_mae}. Run 'sbdd prep'.[/bold red]"
                )
                session.close()
                sys.exit(1)

            plan["ligands"].append(
                {
                    "name": unique_lig_name,
                    "path_sdf": str(lig_sdf),
                    "path_mae": str(lig_mae),
                    "is_native": True,
                    "receptor": r.name,
                    "ligand_title": native_lig.name,
                }
            )

        console.print(
            f"[bold green]✓ 1-to-1 Redocking Mode: Matched {len(active_frames)} receptors directly to their native ligands.[/bold green]"
        )

    else:
        if not ligand_names:
            console.print(
                "[bold red]FATAL: Parameter --ligands is required for cross-docking matrix mode.[/bold red]"
            )
            sys.exit(1)

        ligands = session.query(Ligand).filter(Ligand.name.in_(ligand_names)).all()
        if len(ligands) != len(ligand_names):
            found = {l.name for l in ligands}
            missing = set(ligand_names) - found
            console.print(
                f"[bold red]FATAL: Ligands not found in DB: {missing}[/bold red]"
            )
            sys.exit(1)

        for r in frames:
            sys_obj = session.query(System).filter_by(receptor_frame_id=r.id).first()
            if not sys_obj or not sys_obj.pose:
                console.print(
                    f"[bold red]FATAL: Receptor '{r.name}' has no reference pose! Screening requires a pocket reference.[/bold red]"
                )
                sys.exit(1)

            pose_obj = sys_obj.pose
            native_lig = pose_obj.ligand.name
            pose_p = Path(pose_obj.path).resolve()

            ref_sdf = pose_p.with_suffix(".sdf")
            ref_mae = pose_p.with_suffix(".mae")

            if not ref_sdf.exists() or not ref_mae.exists():
                console.print(
                    f"[bold red]FATAL: Reference pose files missing for '{r.name}'. Need both .sdf and .mae![/bold red]"
                )
                sys.exit(1)

            if not r.path_mae or not r.path_pdb:
                raise ValueError(
                    f"FATAL: ReceptorFrame '{r.name}' does not have both path_mae and path_pdb registered!"
                )

            path_mae = Path(r.path_mae).resolve()
            path_pdb = Path(r.path_pdb).resolve()

            if not path_mae.exists() or not path_pdb.exists():
                raise FileNotFoundError(
                    f"FATAL: Receptor deliverables missing for '{r.name}' on disk"
                )

            if not pose_obj.path_sdf or not pose_obj.path_mae:
                raise ValueError(
                    f"FATAL: Reference Pose '{pose_obj.name}' missing path_sdf or path_mae registration!"
                )

            ref_sdf = Path(pose_obj.path_sdf).resolve()
            ref_mae = Path(pose_obj.path_mae).resolve()

            if not ref_sdf.exists() or not ref_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Reference pose deliverables missing for '{r.name}' on disk"
                )

            plan["receptors"].append(
                {
                    "name": r.name,
                    "path_mae": str(path_mae),
                    "path_pdb": str(path_pdb),
                    "ref_sdf": str(ref_sdf),
                    "ref_mae": str(ref_mae),
                    "native_ligand": native_lig,
                }
            )

        for l in ligands:
            if l.state != "ready" or not l.path_mae or not l.path_sdf:
                console.print(
                    f"[bold red]FATAL: Chemical Ligand '{l.name}' is in state '{l.state}'! Run 'sbdd prep' first.[/bold red]"
                )
                session.close()
                sys.exit(1)

            lig_mae = Path(l.path_mae).resolve()
            lig_sdf = Path(l.path_sdf).resolve()

            if not lig_mae.exists() or not lig_sdf.exists():
                console.print(
                    f"[bold red]FATAL: Prepared files missing for '{l.name}' ({lig_mae}). Run 'sbdd prep' first.[/bold red]"
                )
                session.close()
                sys.exit(1)

            has_crystal_pose = (
                session.query(Pose).filter_by(ligand_id=l.id, origin="crystal").first()
                is not None
            )

            plan["ligands"].append(
                {
                    "name": l.name,
                    "path_sdf": str(lig_sdf),
                    "path_mae": str(lig_mae),
                    "is_native": has_crystal_pose,
                }
            )

        console.print(
            f"[bold cyan]Cross-Docking Matrix Plan generated: {len(frames)} Receptors × {len(ligands)} Ligand Sets[/bold cyan]"
        )

    plan_file = cwd / "dock_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan, f, indent=2)

    # 4. Build Nextflow args
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="screen")
    nf_args.extend(["--dock_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    # 5. Run execution
    ret_code = execute_nextflow(
        "screen", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow docking failed.[/bold red]")
        sys.exit(ret_code)

    # 6. Register results in Database strictly via Plan & Receipt protocol
    num_ligs = len(plan.get("ligands", []))
    default_name = (
        f"redock_{len(frames)}" if is_self else f"screen_{len(frames)}x{num_ligs}"
    )
    poses_name = name if name else default_name

    receipt_file = cwd / "dock_receipt.tsv"
    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not generated by Nextflow. "
            "Docking failed or terminated without producing deliverables."
        )

    res_csv = cwd / "reports" / "ranking.csv"
    if not res_csv.exists() or res_csv.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Ranking summary CSV {res_csv} does not exist. Pipeline failed to summarize docking."
        )

    # Parse receipt: map of (frame_id, engine) -> absolute validated SDF path
    receipt_map = {}
    with open(receipt_file, "r", encoding="utf-8") as rf:
        for line_num, line in enumerate(rf, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                session.close()
                raise ValueError(
                    f"FATAL: Corrupted line {line_num} in {receipt_file}: '{line}'. Expected 3 tab-separated columns."
                )
            r_frame_id, r_engine, r_sdf_path = parts[0], parts[1], parts[2]
            if "/work/" in r_sdf_path or ".nextflow" in r_sdf_path:
                session.close()
                raise ValueError(
                    f"FATAL: Receipt {receipt_file} references scratch path: {r_sdf_path}. "
                    "Must reference published permanent project deliverable."
                )
            sdf_file = Path(r_sdf_path).resolve()
            if not sdf_file.exists() or sdf_file.stat().st_size == 0:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: Published pose file specified in receipt does not exist or is empty: {sdf_file}"
                )
            receipt_map[(r_frame_id, r_engine)] = str(sdf_file)

    import csv
    from rich.table import Table

    table = Table(title="🎯 Redocking & Scoring Summary", show_header=True)
    table.add_column("Receptor", style="cyan")
    table.add_column("Ligand", style="bold")
    table.add_column("GlideScore", justify="right")
    table.add_column("CNN Affinity", justify="right")
    table.add_column("MM-GBSA (kcal/mol)", justify="right")
    table.add_column("Glide RMSD (Å)", justify="center")
    table.add_column("GNINA RMSD (Å)", justify="center")
    table.add_column("Verdict", justify="center")

    poses_count = 0
    with open(res_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("Title", "").strip()
            frame_name = row.get("Receptor", row.get("Frame", "")).strip()
            if not title or not frame_name:
                session.close()
                raise ValueError(
                    f"FATAL: Corrupted ranking row with missing Title or Frame: {row}"
                )

            r_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
            if not r_frame:
                session.close()
                raise ValueError(
                    f"FATAL: ReceptorFrame '{frame_name}' in ranking is not registered in SQLite database!"
                )

            lig = session.query(Ligand).filter_by(name=title).first()
            if not lig:
                session.close()
                raise ValueError(
                    f"FATAL: Ligand '{title}' in ranking is not registered in SQLite database!"
                )

            # Register poses for each executed engine without blind fallbacks
            engines_for_frame = [
                eng for (f_id, eng) in receipt_map.keys() if f_id == frame_name
            ]
            if not engines_for_frame:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: No receipt records found for receptor frame '{frame_name}' in {receipt_file}!"
                )

            for eng in engines_for_frame:
                sdf_path_str = receipt_map[(frame_name, eng)]
                pose_name = f"{frame_name}_{title}_{eng}_dock"
                origin_tag = f"docking_{eng}"

                scores_dict = {}
                if eng == "glide":
                    if row.get("GlideScore") and row["GlideScore"] != "nan":
                        scores_dict["GlideScore"] = float(row["GlideScore"])
                    if (
                        row.get("InSitu_RMSD_Glide_A")
                        and row["InSitu_RMSD_Glide_A"] != "nan"
                    ):
                        scores_dict["RMSD_Glide_A"] = float(row["InSitu_RMSD_Glide_A"])
                elif eng == "gnina":
                    if row.get("Dock_CNNaffinity") and row["Dock_CNNaffinity"] != "nan":
                        scores_dict["CNNaffinity"] = float(row["Dock_CNNaffinity"])
                    if row.get("Dock_CNNscore") and row["Dock_CNNscore"] != "nan":
                        scores_dict["CNNscore"] = float(row["Dock_CNNscore"])
                    if row.get("Dock_VinaAff") and row["Dock_VinaAff"] != "nan":
                        scores_dict["VinaAff"] = float(row["Dock_VinaAff"])
                    if (
                        row.get("InSitu_RMSD_GNINA_A")
                        and row["InSitu_RMSD_GNINA_A"] != "nan"
                    ):
                        scores_dict["RMSD_GNINA_A"] = float(row["InSitu_RMSD_GNINA_A"])

                if row.get("Ultimate_Score") and row["Ultimate_Score"] != "nan":
                    scores_dict["Ultimate_Score"] = float(row["Ultimate_Score"])

                ex_pose = session.query(Pose).filter_by(name=pose_name).first()
                if ex_pose:
                    ex_pose.path = sdf_path_str
                    ex_pose.scores_dict = scores_dict
                    ex_pose.origin = origin_tag
                    ex_pose.state = "ready"
                else:
                    new_p = Pose(
                        ligand_id=lig.id,
                        receptor_frame_id=r_frame.id,
                        name=pose_name,
                        origin=origin_tag,
                        state="ready",
                        path=sdf_path_str,
                    )
                    new_p.scores_dict = scores_dict
                    session.add(new_p)
                poses_count += 1

            # Populate display row for the Rich Summary Table
            gscore = row.get("GlideScore", "-")
            cnn = row.get("Resc_CNNaffinity") or row.get("Dock_CNNaffinity") or "-"
            mmgbsa = row.get("MMGBSA", "-")
            g_rmsd = row.get("InSitu_RMSD_Glide_A", "-")
            n_rmsd = row.get("InSitu_RMSD_GNINA_A", "-")

            disp_glide = _color_rmsd(g_rmsd)
            disp_gnina = _color_rmsd(n_rmsd)

            verdict = "[dim]-[/dim]"
            try:
                vg = float(g_rmsd) if g_rmsd and g_rmsd not in ["nan", "-"] else None
                vn = float(n_rmsd) if n_rmsd and n_rmsd not in ["nan", "-"] else None
                if (vg is not None and vg <= 2.5) or (vn is not None and vn <= 2.5):
                    verdict = "[bold green]✓ SUCCESS[/bold green]"
                elif vg is not None or vn is not None:
                    verdict = "[bold yellow]⚠ DRIFT[/bold yellow]"
            except (ValueError, TypeError):
                pass

            table.add_row(
                frame_name,
                title,
                _format_metric(gscore),
                _format_metric(cnn),
                _format_metric(mmgbsa),
                disp_glide,
                disp_gnina,
                verdict,
            )

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Docking verified via receipt! Registered {poses_count} pose(s) in DB as '{poses_name}'.[/bold green]\n"
    )
    console.print(table)
