"""
File: cli/commands/dock.py
Description: CLI handler for 'sbdd dock'. Reads DB to construct Cross-Docking plan.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Receptor, Ligand, PoseSet, System, Pose
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_dock(args):
    cwd = Path.cwd()
    session = init_db()
    is_self = getattr(args, "self_dock", False)

    # 1. Fetch requested Receptors
    if getattr(args, "receptors", None):
        rec_names = [n.strip() for n in args.receptors.split(",")]
        receptors = session.query(Receptor).filter(Receptor.name.in_(rec_names)).all()
        if len(receptors) != len(rec_names):
            found = {r.name for r in receptors}
            missing = set(rec_names) - found
            console.print(f"[bold red]FATAL: Receptors not found in DB: {missing}[/bold red]")
            sys.exit(1)
    else:
        if is_self:
            receptors = session.query(Receptor).filter_by(status="ready").all()
        else:
            console.print("[bold red]FATAL: Please specify --receptors or use --self / --redock mode.[/bold red]")
            sys.exit(1)

    for r in receptors:
        if r.status != "ready":
            console.print(f"[bold red]FATAL: Receptor '{r.name}' is not ready! Run 'sbdd prep' first.[/bold red]")
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
        for r in receptors:
            # Walk up to root target to resolve the canonical native ligand
            curr = r
            while curr.parent_name:
                curr = session.query(Receptor).filter_by(name=curr.parent_name).first()
            root_target = curr

            root_sys = session.query(System).filter_by(receptor_name=root_target.name).first()
            if not root_sys or not root_sys.pose or not root_sys.pose.ligand:
                console.print(f"[bold red]FATAL: No native ligand found in database for receptor '{r.name}' (Root Target: {root_target.name}).[/bold red]")
                sys.exit(1)

            native_lig = root_sys.pose.ligand
            ligand_set.add(native_lig)

            # Reference pose for grid positioning & RMSD
            ref_pose_path = "NO_FILE"
            r_sys = session.query(System).filter_by(receptor_name=r.name).first()
            if r_sys and r_sys.pose and Path(r_sys.pose.path).exists():
                pose_p = Path(r_sys.pose.path)
                if pose_p.suffix == ".sdf" and pose_p.with_suffix(".mae").exists():
                    ref_pose_path = str(pose_p.with_suffix(".mae"))
                else:
                    ref_pose_path = str(pose_p)

            plan["receptors"].append({
                "name": r.name,
                "path": r.path,
                "ref_pose_path": ref_pose_path,
                "native_ligand": root_target.meta_dict.get("native_ligand", ""),
            })

            plan["pairs"].append({
                "receptor": r.name,
                "ligand": native_lig.name
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
        ex_pose = session.query(PoseSet).filter_by(name=poses_name).first()
        if ex_pose:
            ex_pose.path_csv = str(res_csv)
        else:
            new_pose = PoseSet(
                name=poses_name,
                path_sdf="multiple_files",
                path_csv=str(res_csv),
                status="ready",
            )
            session.add(new_pose)
            
        # Bulk insert top poses if single receptor (Cross-dock bulk mapping is deferred)
        # Bulk insert all docking poses into DB across all receptors
        poses_to_insert = []
        import csv
        try:
            with open(res_csv, 'r', encoding='utf-8', errors='ignore') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    title = row.get("Title", "").strip()
                    rec_name = row.get("Receptor", row.get("Frame", "")).strip()
                    if not title or not rec_name:
                        continue

                    pose_name = f"{title}_{rec_name}_dock"
                    sdf_path = cwd / "docking" / rec_name / f"{rec_name}_dock.sdf"
                    if not sdf_path.exists():
                        sdf_path = cwd / "docking" / f"{rec_name}_dock.sdf"

                    scores_dict = {
                        k: float(v) for k, v in row.items()
                        if k not in ["Rank", "Title", "Receptor", "Frame"] and v and v != "nan" and v != ""
                    }

                    ex = session.query(Pose).filter_by(name=pose_name).first()
                    if not ex:
                        lig = session.query(Ligand).filter((Ligand.name == title) | (Ligand.name.like(f"%{title}%"))).first()
                        lig_db_name = lig.name if lig else title

                        poses_to_insert.append(Pose(
                            name=pose_name,
                            ligand_name=lig_db_name,
                            receptor_name=rec_name,
                            origin="docking",
                            path=str(sdf_path),
                            status="ready",
                            scores=json.dumps(scores_dict)
                        ))
            if poses_to_insert:
                session.bulk_save_objects(poses_to_insert)
        except Exception as e:
            console.print(f"[yellow]Warning: Could not bulk insert poses: {e}[/yellow]")

        session.commit()
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
