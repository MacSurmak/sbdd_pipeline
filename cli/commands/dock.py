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
            console.print(f"[bold red]FATAL: Target '{args.target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif getattr(args, "frames", None) or getattr(args, "receptors", None):
        req = getattr(args, "frames", None) or getattr(args, "receptors", None)
        f_names = [f.strip() for f in req.split(",") if f.strip()]
        query = query.filter(ReceptorFrame.name.in_(f_names))

    frames = query.all()
    if not frames:
        console.print("[bold red]FATAL: No ready ReceptorFrames found to dock. Run 'sbdd prep' first![/bold red]")
        sys.exit(1)

    # 2. Build Plan (Self 1-to-1 vs Matrix Screening)
    plan = {
        "mode": "self" if is_self else "matrix",
        "receptors": [],
        "ligands": [],
        "pairs": []
    }

    if is_self:
        ligand_set = set()
        for r in frames:
            # Resolve native ligand directly from target's systems
            target_sys = session.query(System).filter(System.target_id == r.target_id, System.pose_id.isnot(None)).first()
            if not target_sys or not target_sys.pose:
                console.print(f"[bold red]FATAL: No reference pose found for frame '{r.name}' in target '{r.target.name}'.[/bold red]")
                sys.exit(1)

            native_pose = target_sys.pose
            native_lig = native_pose.ligand
            ligand_set.add((native_lig, native_pose.path))

            ref_pose_path = native_pose.path
            pose_p = Path(ref_pose_path)
            if pose_p.suffix == ".sdf" and pose_p.with_suffix(".mae").exists():
                ref_pose_path = str(pose_p.with_suffix(".mae"))

            plan["receptors"].append({
                "name": r.name,
                "path": r.path,
                "ref_pose_path": ref_pose_path,
                "native_ligand": native_lig.name,
            })

            plan["pairs"].append({
                "receptor": r.name,
                "ligand": native_lig.name
            })

        for l_obj, l_path in sorted(list(ligand_set), key=lambda x: x[0].name):
            plan["ligands"].append({
                "name": l_obj.name,
                "path": str(Path(l_path).resolve()),
                "status": "ready",
                "is_native": True
            })

        for l in sorted(list(ligand_set), key=lambda x: x.name):
            plan["ligands"].append({
                "name": l.name,
                "path": str(Path(l.path).resolve()),
                "status": l.status,
                "is_native": True
            })

        console.print(f"[bold green]✓ 1-to-1 Redocking Mode: Matched {len(plan['pairs'])} receptors directly to their native ligands (0 cross-docking jobs).[/bold green]")

    else:
        if not getattr(args, "ligands", None):
            console.print("[bold red]FATAL: Parameter --ligands is required for cross-docking matrix mode.[/bold red]")
            sys.exit(1)

        lig_names = [n.strip() for n in args.ligands.split(",")]
        ligands = session.query(Ligand).filter(Ligand.name.in_(lig_names)).all()
        if len(ligands) != len(lig_names):
            found = {l.name for l in ligands}
            missing = set(lig_names) - found
            console.print(f"[bold red]FATAL: Ligands not found in DB: {missing}[/bold red]")
            sys.exit(1)

        for r in receptors:
            ref_pose_path = "NO_FILE"
            sys_obj = session.query(System).filter_by(receptor_name=r.name).first()
            if sys_obj and sys_obj.pose_name:
                pose_obj = session.query(Pose).filter_by(name=sys_obj.pose_name).first()
                if pose_obj and Path(pose_obj.path).exists():
                    pose_p = Path(pose_obj.path)
                    if pose_p.suffix == ".sdf" and pose_p.with_suffix(".mae").exists():
                        ref_pose_path = str(pose_p.with_suffix(".mae"))
                    else:
                        ref_pose_path = str(pose_p)

            plan["receptors"].append({
                "name": r.name,
                "path": r.path,
                "ref_pose_path": ref_pose_path,
                "native_ligand": r.meta_dict.get("native_ligand", ""),
            })

        for l in ligands:
            has_crystal_pose = (
                session.query(Pose).filter_by(ligand_name=l.name, origin="crystal").first() is not None
            )
            plan["ligands"].append({
                "name": l.name,
                "path": str(Path(l.path).resolve()),
                "status": l.status,
                "is_native": has_crystal_pose,
            })

        console.print(f"[bold cyan]Cross-Docking Matrix Plan generated: {len(receptors)} Receptors × {len(ligands)} Ligand Sets[/bold cyan]")

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
    default_name = f"redock_{len(receptors)}" if is_self else f"screen_{len(receptors)}x{num_ligs}"
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

                    r_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
                    if not r_frame:
                        continue

                    lig = session.query(Ligand).filter_by(name=title).first()
                    if not lig:
                        lig = session.query(Ligand).filter(Ligand.name.like(f"%{title}%")).first()
                    if not lig:
                        continue

                    pose_name = f"{frame_name}_{title}_dock"
                    sdf_path = cwd / "docking" / frame_name / f"{frame_name}_dock.sdf"
                    if not sdf_path.exists():
                        sdf_path = cwd / "docking" / f"{frame_name}_dock.sdf"

                    scores_dict = {
                        k: float(v) for k, v in row.items()
                        if k not in ["Rank", "Title", "Receptor", "Frame"] and v and v != "nan" and v != ""
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
            console.print(f"[yellow]Warning: Could not record docking poses in DB: {e}[/yellow]")
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

            with open(res_csv, 'r', encoding='utf-8', errors='ignore') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rec = row.get("Receptor", row.get("Frame", "-"))
                    lig = row.get("Title", "-")
                    gscore = row.get("GlideScore", "-")
                    cnn = row.get("Resc_CNNaffinity") or row.get("Dock_CNNaffinity") or "-"
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

                    table.add_row(rec, lig, _fmt(gscore), _fmt(cnn), _fmt(mmgbsa), rmsd_display, verdict)
            console.print(table)
        except Exception as e:
            console.print(f"[yellow]Warning: Could not render console table: {e}[/yellow]")

    session.close()
