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


def cmd_dock(args):
    cwd = Path.cwd()
    session = init_db()
    is_self = getattr(args, "self_dock", False)

    # 1. Fetch requested Ready Frames
    query = session.query(ReceptorFrame).filter_by(state="ready")
    if getattr(args, "target", None):
        target_obj = session.query(Target).filter_by(name=args.target).first()
        if not target_obj:
            console.print(
                f"[bold red]FATAL: Target '{args.target}' not found.[/bold red]"
            )
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif getattr(args, "frames", None) or getattr(args, "receptors", None):
        req = getattr(args, "frames", None) or getattr(args, "receptors", None)
        f_names = [f.strip() for f in req.split(",") if f.strip()]
        query = query.filter(ReceptorFrame.name.in_(f_names))

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
            # Strictly resolve existing MAE and PDB files for the receptor
            frame_p = Path(r.path).resolve()
            path_mae = frame_p.with_suffix(".mae")
            path_pdb = frame_p.with_suffix(".pdb")

            if not path_mae.exists():
                console.print(
                    f"[bold red]FATAL: Receptor MAE not found for frame '{r.name}': {path_mae}[/bold red]"
                )
                sys.exit(1)
            if not path_pdb.exists():
                console.print(
                    f"[bold red]FATAL: Receptor PDB not found for frame '{r.name}': {path_pdb}[/bold red]"
                )
                sys.exit(1)

            # Strictly resolve existing SDF and MAE files for reference pose
            pose_p = Path(native_pose.path).resolve()
            pose_sdf = pose_p.with_suffix(".sdf")
            pose_mae = pose_p.with_suffix(".mae")

            if not pose_sdf.exists():
                console.print(
                    f"[bold red]FATAL: Reference pose SDF not found for frame '{r.name}': {pose_sdf}[/bold red]"
                )
                sys.exit(1)
            if not pose_mae.exists():
                console.print(
                    f"[bold red]FATAL: Reference pose MAE not found for frame '{r.name}': {pose_mae}[/bold red]"
                )
                sys.exit(1)

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
            plan["ligands"].append(
                {
                    "name": unique_lig_name,
                    "path": str(Path(native_pose.path).resolve()),
                    "status": "ready" if native_lig.batch_id == "crystal" else "unprepared",
                    "is_native": True,
                    "receptor": r.name,
                    "ligand_title": native_lig.name,
                }
            )

        console.print(
            f"[bold green]✓ 1-to-1 Redocking Mode: Matched {len(active_frames)} receptors directly to their native ligands.[/bold green]"
        )

    else:
        if not getattr(args, "ligands", None):
            console.print(
                "[bold red]FATAL: Parameter --ligands is required for cross-docking matrix mode.[/bold red]"
            )
            sys.exit(1)

        lig_names = [n.strip() for n in args.ligands.split(",")]
        ligands = session.query(Ligand).filter(Ligand.name.in_(lig_names)).all()
        if len(ligands) != len(lig_names):
            found = {l.name for l in ligands}
            missing = set(lig_names) - found
            console.print(
                f"[bold red]FATAL: Ligands not found in DB: {missing}[/bold red]"
            )
            sys.exit(1)

        for r in frames:
            sys_obj = session.query(System).filter_by(receptor_frame_id=r.id).first()
            if not sys_obj or not sys_obj.pose:
                console.print(f"[bold red]FATAL: Receptor '{r.name}' has no reference pose! Screening requires a pocket reference.[/bold red]")
                sys.exit(1)
                
            pose_obj = sys_obj.pose
            native_lig = pose_obj.ligand.name
            pose_p = Path(pose_obj.path).resolve()
            
            ref_sdf = pose_p.with_suffix(".sdf")
            ref_mae = pose_p.with_suffix(".mae")
            
            if not ref_sdf.exists() or not ref_mae.exists():
                console.print(f"[bold red]FATAL: Reference pose files missing for '{r.name}'. Need both .sdf and .mae![/bold red]")
                sys.exit(1)

            frame_p = Path(r.path).resolve()
            path_mae = frame_p.with_suffix(".mae")
            path_pdb = frame_p.with_suffix(".pdb")
            
            if not path_mae.exists() or not path_pdb.exists():
                console.print(f"[bold red]FATAL: Receptor files missing for '{r.name}'. Need both .pdb and .mae![/bold red]")
                sys.exit(1)

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
            has_crystal_pose = (
                session.query(Pose)
                .filter_by(ligand_id=l.id, origin="crystal")
                .first()
                is not None
            )
            # Use 'ready' status for standalone chemical ligands that don't need prep
            plan["ligands"].append(
                {
                    "name": l.name,
                    "path": str(Path(l.path).resolve()),
                    "status": "ready" if l.batch_id == "crystal" else "unprepared",
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

    if getattr(args, "engine", None):
        if "--dock_engine" in nf_args:
            idx = nf_args.index("--dock_engine")
            nf_args[idx + 1] = args.engine
        else:
            nf_args.extend(["--dock_engine", args.engine])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    # 5. Run execution
    ret_code = execute_nextflow(
        "screen", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow docking failed.[/bold red]")
        sys.exit(ret_code)

    # 6. Register results in Database
    num_ligs = len(plan.get("ligands", []))
    default_name = (
        f"redock_{len(frames)}" if is_self else f"screen_{len(frames)}x{num_ligs}"
    )
    poses_name = args.name if getattr(args, "name", None) else default_name
    res_csv = cwd / "reports" / "ranking.csv"
    res_sdf = (
        cwd / "docking" / f"{poses_name}_poses.sdf"
    )  # We will implement SDF merging later if needed, tracking CSV for now

    if res_csv.exists():
        import csv

        poses_count = 0
        try:
            with open(res_csv, "r", encoding="utf-8", errors="ignore") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    title = row.get("Title", "").strip()
                    frame_name = row.get("Receptor", row.get("Frame", "")).strip()
                    if not title or not frame_name:
                        continue

                    r_frame = (
                        session.query(ReceptorFrame).filter_by(name=frame_name).first()
                    )
                    if not r_frame:
                        continue

                    lig = session.query(Ligand).filter_by(name=title).first()
                    if not lig:
                        lig = (
                            session.query(Ligand)
                            .filter(Ligand.name.like(f"%{title}%"))
                            .first()
                        )
                    if not lig:
                        continue

                    pose_name = f"{frame_name}_{title}_dock"
                    
                    sdf_path = None
                    for cand in [
                        cwd / "docking" / frame_name / f"{frame_name}_dock.sdf",
                        cwd / "docking" / frame_name / f"{frame_name}_gnina_docked_all.sdf",
                        cwd / "docking" / f"{frame_name}_dock.sdf",
                        cwd / "docking" / f"{frame_name}_gnina_docked_all.sdf",
                    ]:
                        if cand.exists():
                            sdf_path = cand
                            break

                    if not sdf_path:
                        console.print(f"[bold red]FATAL: Could not find docked SDF for {pose_name}. Pipeline failed to produce poses![/bold red]")
                        sys.exit(1)

                    scores_dict = {
                        k: float(v)
                        for k, v in row.items()
                        if k not in ["Rank", "Title", "Receptor", "Frame"]
                        and v
                        and v != "nan"
                        and v != ""
                    }

                    ex = session.query(Pose).filter_by(name=pose_name).first()
                    if ex:
                        ex.path = str(sdf_path.resolve())
                        ex.scores_dict = scores_dict
                        ex.state = "ready"
                    else:
                        new_p = Pose(
                            ligand_id=lig.id,
                            receptor_frame_id=r_frame.id,
                            name=pose_name,
                            origin="docking",
                            state="ready",
                            path=str(sdf_path.resolve()),
                        )
                        new_p.scores_dict = scores_dict
                        session.add(new_p)
                    poses_count += 1
            session.commit()
        except Exception as e:
            console.print(
                f"[yellow]Warning: Could not record docking poses in DB: {e}[/yellow]"
            )
        console.print(
            f"\n[bold green]✓ Docking complete! Results registered as PoseSet '{poses_name}'.[/bold green]\n"
        )

        # Render Rich Summary Table
        try:
            from rich.table import Table

            table = Table(title="🎯 Redocking & Scoring Summary", show_header=True)
            table.add_column("Receptor", style="cyan")
            table.add_column("Ligand", style="bold")
            table.add_column("GlideScore", justify="right")
            table.add_column("CNN Affinity", justify="right")
            table.add_column("MM-GBSA (kcal/mol)", justify="right")
            table.add_column("InSitu RMSD (Å)", justify="center")
            table.add_column("Verdict", justify="center")

            with open(res_csv, "r", encoding="utf-8", errors="ignore") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rec = row.get("Receptor", row.get("Frame", "-"))
                    lig = row.get("Title", "-")
                    gscore = row.get("GlideScore", "-")
                    cnn = (
                        row.get("Resc_CNNaffinity")
                        or row.get("Dock_CNNaffinity")
                        or "-"
                    )
                    mmgbsa = row.get("MMGBSA", "-")
                    rmsd_str = row.get("InSitu_RMSD_A", "-")

                    verdict = "[dim]-[/dim]"
                    rmsd_display = "-"
                    if rmsd_str and rmsd_str != "nan" and rmsd_str != "-":
                        try:
                            val = float(rmsd_str)
                            if val <= 2.0:
                                rmsd_display = f"[bold green]{val:.2f}[/bold green]"
                                verdict = "[bold green]✓ SUCCESS[/bold green]"
                            else:
                                rmsd_display = f"[bold red]{val:.2f}[/bold red]"
                                verdict = "[bold yellow]⚠ DRIFT[/bold yellow]"
                        except ValueError:
                            rmsd_display = rmsd_str

                    def _fmt(val):
                        try:
                            return f"{float(val):.2f}"
                        except (ValueError, TypeError):
                            return str(val) if val else "-"

                    table.add_row(
                        rec,
                        lig,
                        _fmt(gscore),
                        _fmt(cnn),
                        _fmt(mmgbsa),
                        rmsd_display,
                        verdict,
                    )
            console.print(table)
        except Exception as e:
            console.print(
                f"[yellow]Warning: Could not render console table: {e}[/yellow]"
            )

    session.close()
