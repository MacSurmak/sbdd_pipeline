"""
File: cli/commands/traj.py
Description: Clean 3NF CLI handlers for trajectory solvent stripping, metrics, and clustering.
Enforces stream cleaning in worker scratch without raw file duplication on disk.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from rich.table import Table
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System, Trajectory
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def _query_trajectories(session, args, target_status: str = "raw") -> list[Trajectory]:
    """
    Queries trajectories needing action based on DB status.
    By default queries all trajectories with status == target_status across all systems.
    Allows optional scope restriction via --systems, --target, or --traj.
    """
    query = session.query(Trajectory)

    if getattr(args, "systems", None):
        sys_names = [s.strip() for s in args.systems.split(",") if s.strip()]
        systems = session.query(System).filter(System.name.in_(sys_names)).all()
        sys_ids = [s.id for s in systems]
        query = query.filter(Trajectory.system_id.in_(sys_ids))
    elif getattr(args, "target", None):
        target_obj = session.query(Target).filter_by(name=args.target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{args.target}' not found.[/bold red]")
            sys.exit(1)
        systems = session.query(System).filter_by(target_id=target_obj.id).all()
        sys_ids = [s.id for s in systems]
        query = query.filter(Trajectory.system_id.in_(sys_ids))
    elif getattr(args, "traj", None):
        traj_names = [t.strip() for t in args.traj.split(",") if t.strip()]
        query = query.filter(Trajectory.name.in_(traj_names))

    # Strict state filter: process only trajectories that actually need it
    if target_status:
        query = query.filter(Trajectory.status == target_status)

    trajs = query.all()
    return trajs


def cmd_clean_traj(args):
    """
    Cleans raw trajectories (PBC correction, fit, solvent strip).
    Streams chunks into worker scratch without duplicating raw files on storage.
    Automatically skips already cleaned trajectories.
    """
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(session, args, target_status="raw")

    if not trajs:
        console.print("[bold green]✓ All trajectories are already cleaned! Nothing to do.[/bold green]")
        session.close()
        return

    plan_list = []
    for t in trajs:
        raw_dir = Path(t.raw_dir).resolve()
        chunks = sorted(
            raw_dir.glob("md_*.xtc"),
            key=lambda p: int("".join(filter(str.isdigit, p.stem)) or 0)
        )
        if not chunks:
            # Fallback if raw directory contains single xtc
            chunks = list(raw_dir.glob("*.xtc"))

        if not chunks:
            console.print(f"[bold red]FATAL: No raw trajectory chunks found in {raw_dir}[/bold red]")
            sys.exit(1)

        clean_dir = raw_dir.parent / "clean"
        plan_list.append(
            {
                "name": t.name,
                "rec_id": t.system.name,
                "rep": t.replicate,
                "xtcs": [str(c.resolve()) for c in chunks],
                "tpr": str(Path(t.path_tpr).resolve()),
                "clean_dir": str(clean_dir.resolve()),
            }
        )

    plan_file = cwd / "clean_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(f"[bold cyan]Cleaning {len(trajs)} trajectory(ies) (Stream merge in worker scratch -> Fit -> Strip)...[/bold cyan]")

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--clean_plan", str(plan_file.resolve())])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow trajectory cleaning failed.[/bold red]")
        sys.exit(ret_code)

    # Update database records
    for item in plan_list:
        clean_dir = Path(item["clean_dir"])
        clean_xtc = clean_dir / "trajectory_clean.xtc"
        ref_pdb = clean_dir / "complex_ref.pdb"

        if clean_xtc.exists() and ref_pdb.exists():
            t_obj = session.query(Trajectory).filter_by(name=item["name"]).first()
            if t_obj:
                t_obj.path_clean_xtc = str(clean_xtc.resolve())
                t_obj.path_ref_pdb = str(ref_pdb.resolve())
                t_obj.status = "cleaned"

    session.commit()
    session.close()

    console.print(f"\n[bold green]✓ Successfully cleaned {len(trajs)} trajectory(ies) without raw file duplication![/bold green]")
    console.print("[dim]Clean trajectories saved to clean/ (130 MB). Chunks in raw/ preserved for restarts.[/dim]")


def cmd_metrics(args):
    """Calculates RMSD (Protein, Ligand, and Cofactors), RMSF, and Contacts on cleaned trajectories."""
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(session, args, target_status="cleaned")

    if not trajs:
        console.print("[yellow]No cleaned trajectories found to analyze. Run 'sbdd traj-clean' first![/yellow]")
        session.close()
        return

    plan_list = []
    for t in trajs:
        if not t.path_clean_xtc or not t.path_ref_pdb:
            continue

        sys_obj = t.system
        lig_code = sys_obj.pose.ligand.name if sys_obj.pose else ""
        cofactors = sys_obj.receptor_frame.cofactors_list

        clean_path = Path(t.path_clean_xtc)
        metrics_dir = clean_path.parent.parent / "metrics"

        plan_list.append(
            {
                "name": t.name,
                "rec_id": sys_obj.name,
                "rep": t.replicate,
                "clean_xtc": str(clean_path.resolve()),
                "ref_pdb": str(Path(t.path_ref_pdb).resolve()),
                "ligand": lig_code,
                "cofactors": cofactors,
                "metrics_dir": str(metrics_dir.resolve()),
            }
        )

    plan_file = cwd / "metrics_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(f"[bold cyan]Computing metrics for {len(trajs)} clean trajectory(ies)...[/bold cyan]")

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--metrics_plan", str(plan_file.resolve())])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow metrics calculation failed.[/bold red]")
        sys.exit(ret_code)

    # Render results table
    table = Table(title="Molecular Dynamics Trajectory Metrics", show_header=True)
    table.add_column("Trajectory", style="cyan")
    table.add_column("Prot BB RMSD (Å)", justify="center")
    table.add_column("Ligand RMSD (Å)", justify="center", style="green")
    table.add_column("Cofactors RMSD (Å)", justify="center", style="yellow")
    table.add_column("Plot Deliverable", style="dim")

    for item in plan_list:
        summary_f = Path(item["metrics_dir"]) / "md_metrics_summary.json"
        plot_f = Path(item["metrics_dir"]) / "rmsd_plot.png"
        if summary_f.exists():
            data = json.loads(summary_f.read_text())
            p_mean = data.get("Prot_BB_RMSD_Mean_A", "-")
            p_std = data.get("Prot_BB_RMSD_Std_A", "-")
            prot_str = f"{p_mean} ± {p_std}" if p_mean != "-" else "-"

            lig_str = "-"
            cof_strs = []
            for m_key, m_val in data.get("Molecules", {}).items():
                if m_key.startswith("Lig_"):
                    lig_str = f"{m_val.get('Mean_A', '-')} ± {m_val.get('Std_A', '-')}"
                elif m_key.startswith("Cofactor_"):
                    c_name = m_key.replace("Cofactor_", "").replace("_RMSD_A", "")
                    cof_strs.append(f"{c_name}: {m_val.get('Mean_A', '-')}")

            cof_display = ", ".join(cof_strs) if cof_strs else "None"
            rel_plot = str(plot_f.relative_to(cwd)) if plot_f.exists() else "-"
            table.add_row(item["name"], prot_str, lig_str, cof_display, rel_plot)

    console.print(table)
    session.close()


def cmd_cluster(args):
    """
    Clusters clean trajectory ensembles and registers each centroid as a
    peer ReceptorFrame (state='ready') in the parent Target.
    Automatically skips systems whose centroids are already generated.
    """
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(session, args, target_status="cleaned")

    if not trajs:
        console.print("[yellow]No cleaned trajectories found to cluster. Run 'sbdd traj-clean' first![/yellow]")
        session.close()
        return

    for t in trajs:
        if t.status != "cleaned" or not t.path_clean_xtc or not t.path_ref_pdb:
            console.print(f"[bold red]FATAL: Trajectory '{t.name}' must be cleaned first. Run 'sbdd traj-clean'![/bold red]")
            sys.exit(1)

    # Group trajectories by System
    by_sys = {}
    for t in trajs:
        by_sys.setdefault(t.system_id, []).append(t)

    plan_list = []
    sys_map = {}

    for sys_id, t_list in by_sys.items():
        sys_obj = session.query(System).get(sys_id)
        task_id = f"{sys_obj.name}_cluster"
        sys_map[task_id] = sys_obj

        clean_xtcs = [t.path_clean_xtc for t in t_list]
        ref_pdb = t_list[0].path_ref_pdb
        lig_code = sys_obj.pose.ligand.name if sys_obj.pose else ""

        plan_list.append(
            {
                "task_id": task_id,
                "system": sys_obj.name,
                "clean_xtcs": clean_xtcs,
                "ref_pdb": ref_pdb,
                "ligand": lig_code,
                "cutoff": args.cutoff,
                "method": args.method,
                "dt": args.dt,
                "group": args.group,
                "max_clusters": getattr(args, "max_clusters", 5),
            }
        )

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(f"[bold cyan]Clustering trajectories across {len(by_sys)} System(s)...[/bold cyan]")

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--cluster_plan", str(plan_file.resolve())])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow clustering failed.[/bold red]")
        sys.exit(ret_code)

    receipt_file = cwd / "cluster_receipt.tsv"
    if not receipt_file.exists():
        console.print("[bold red]FATAL: cluster_receipt.tsv not found![/bold red]")
        sys.exit(1)

    registered_frames = 0
    for line in receipt_file.read_text().strip().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        cluster_id = parts[0]
        pdb_path = parts[1]
        pose_path = parts[2] if len(parts) > 2 else "NO_POSE"

        if "_pose" in cluster_id:
            continue

        parent_sys_key = cluster_id.split("_centroid_")[0]
        parent_sys = sys_map.get(parent_sys_key)
        if not parent_sys:
            continue

        # A. Register Centroid as a peer ReceptorFrame in Target
        centroid_frame_name = f"{parent_sys.name}_{cluster_id.split('_')[-2]}_{cluster_id.split('_')[-1]}"
        ex_frame = session.query(ReceptorFrame).filter_by(name=centroid_frame_name).first()
        if not ex_frame:
            new_frame = ReceptorFrame(
                target_id=parent_sys.target_id,
                name=centroid_frame_name,
                source_type="md_centroid",
                state="ready",
                parent_frame_id=parent_sys.receptor_frame_id,
                path=str(Path(pdb_path).resolve()),
            )
            new_frame.cofactors_list = parent_sys.receptor_frame.cofactors_list
            session.add(new_frame)
            session.flush()
            frame_obj = new_frame
        else:
            frame_obj = ex_frame

        # B. Register in-situ Pose if extracted from cluster
        pose_obj = None
        if pose_path != "NO_POSE" and Path(pose_path).exists() and parent_sys.pose:
            cluster_pose_name = f"{centroid_frame_name}_pose"
            ex_pose = session.query(Pose).filter_by(name=cluster_pose_name).first()
            if not ex_pose:
                new_pose = Pose(
                    ligand_id=parent_sys.pose.ligand_id,
                    receptor_frame_id=frame_obj.id,
                    name=cluster_pose_name,
                    origin="md_cluster",
                    state="ready",
                    path=str(Path(pose_path).resolve()),
                )
                session.add(new_pose)
                session.flush()
                pose_obj = new_pose
            else:
                pose_obj = ex_pose

        # C. Register Ready System for this Centroid
        cluster_sys_name = f"{centroid_frame_name}_holo" if pose_obj else f"{centroid_frame_name}_apo"
        ex_sys = session.query(System).filter_by(name=cluster_sys_name).first()
        if not ex_sys:
            new_sys = System(
                target_id=parent_sys.target_id,
                receptor_frame_id=frame_obj.id,
                pose_id=pose_obj.id if pose_obj else None,
                name=cluster_sys_name,
                state="ready",
                origin="md_cluster",
            )
            session.add(new_sys)
        registered_frames += 1

    session.commit()
    session.close()

    console.print(f"\n[bold green]✓ Clustering complete! Registered {registered_frames} centroid frame(s) as ready ReceptorFrames.[/bold green]")
    console.print("[dim]These frames can now be used directly for docking. View with 'sbdd tree'.[/dim]")
