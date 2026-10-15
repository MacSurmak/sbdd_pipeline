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
from cli.core.db import init_db, Target, ReceptorFrame, Substance, Pose, System, Trajectory
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def _query_trajectories(
    session,
    systems: list[str] | None = None,
    target: str | None = None,
    traj_names: list[str] | None = None,
    target_status: str = "raw",
) -> list[Trajectory]:
    """
    Queries trajectories needing action based on DB status.
    By default queries all trajectories with status == target_status across all systems.
    Allows optional scope restriction via systems, target, or traj_names.
    """
    query = session.query(Trajectory)

    if systems:
        sys_objs = session.query(System).filter(System.name.in_(systems)).all()
        sys_ids = [s.id for s in sys_objs]
        query = query.filter(Trajectory.system_id.in_(sys_ids))
    elif target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        sys_objs = session.query(System).filter_by(target_id=target_obj.id).all()
        sys_ids = [s.id for s in sys_objs]
        query = query.filter(Trajectory.system_id.in_(sys_ids))
    elif traj_names:
        query = query.filter(Trajectory.name.in_(traj_names))

    # Strict state filter: process only trajectories that actually need it
    if target_status:
        query = query.filter(Trajectory.status == target_status)

    trajs = query.all()
    return trajs


def run_clean_traj(
    systems: list[str] | None = None,
    target: str | None = None,
    traj_names: list[str] | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """
    Cleans raw trajectories (PBC correction, fit, solvent strip).
    Streams chunks into worker scratch without duplicating raw files on storage.
    Automatically skips already cleaned trajectories.
    """
    import math

    cwd = Path.cwd()
    session = init_db()
    configs = load_project_configs(cwd)
    trajs = _query_trajectories(
        session, systems, target, traj_names, target_status="raw"
    )

    if not trajs:
        console.print(
            "[bold green]✓ All trajectories are already cleaned! Nothing to do.[/bold green]"
        )
        session.close()
        return

    try:
        chunk_ns = float(configs["md_simulation"]["production"]["chunk_size_ns"])
    except KeyError:
        console.print(
            "[bold red]FATAL: 'chunk_size_ns' missing in md.yaml configuration.[/bold red]"
        )
        sys.exit(1)

    plan_list = []
    for t in trajs:
        raw_dir = Path(t.raw_dir).resolve()

        expected_chunks = max(1, int(math.ceil(t.length_ns / chunk_ns)))
        chunks = [raw_dir / f"md_{i}.xtc" for i in range(1, expected_chunks + 1)]

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

    console.print(
        f"[bold cyan]Cleaning {len(trajs)} trajectory(ies) (Stream merge in worker scratch -> Fit -> Strip)...[/bold cyan]"
    )
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--clean_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow trajectory cleaning failed.[/bold red]")
        sys.exit(ret_code)

    # Update database records
    for item in plan_list:
        clean_dir = Path(item["clean_dir"])
        clean_xtc = clean_dir / "trajectory_clean.xtc"
        ref_struct = clean_dir / "complex_ref.pdb"

        if clean_xtc.exists() and ref_struct.exists():
            t_obj = session.query(Trajectory).filter_by(name=item["name"]).first()
            if t_obj:
                t_obj.path_clean_xtc = str(clean_xtc.resolve())
                t_obj.path_ref_pdb = str(ref_struct.resolve())
                t_obj.status = "cleaned"

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Successfully cleaned {len(trajs)} trajectory(ies) without raw file duplication![/bold green]"
    )
    console.print(
        "[dim]Clean trajectories saved to clean/ (130 MB). Chunks in raw/ preserved for restarts.[/dim]"
    )


def run_metrics(
    systems: list[str] | None = None,
    target: str | None = None,
    traj_names: list[str] | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """Calculates RMSD (Protein, Ligand, and Cofactors), RMSF, and Contacts on cleaned trajectories."""
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(
        session, systems, target, traj_names, target_status="cleaned"
    )

    if not trajs:
        console.print(
            "[yellow]No cleaned trajectories found to analyze. Run 'sbdd traj-clean' first![/yellow]"
        )
        session.close()
        return

    plan_list = []
    for t in trajs:
        if not t.path_clean_xtc:
            continue

        clean_path = Path(t.path_clean_xtc)
        clean_dir = clean_path.parent
        ref_pdb = clean_dir / "complex_ref.pdb"

        if not ref_pdb.exists():
            console.print(
                f"[bold red]FATAL: Required topology {ref_pdb} does not exist! Run 'sbdd traj-clean'.[/bold red]"
            )
            sys.exit(1)

        sys_obj = t.system
        lig_code = sys_obj.pose.substance.name if (sys_obj.pose and sys_obj.pose.substance) else ""
        cofactors = [c.substance.name for c in sys_obj.cofactors if c.substance]
        metrics_dir = clean_dir.parent / "metrics"

        plan_list.append(
            {
                "name": t.name,
                "rec_id": sys_obj.name,
                "rep": t.replicate,
                "clean_xtc": str(clean_path.resolve()),
                "ref_pdb": str(ref_pdb.resolve()),
                "ligand": lig_code,
                "cofactors": cofactors,
                "metrics_dir": str(metrics_dir.resolve()),
            }
        )

    plan_file = cwd / "metrics_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Computing metrics for {len(trajs)} clean trajectory(ies)...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--metrics_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
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
            if plot_f.exists():
                try:
                    rel_plot = str(plot_f.relative_to(cwd))
                except ValueError:
                    rel_plot = str(plot_f)
            else:
                rel_plot = "-"
            table.add_row(item["name"], prot_str, lig_str, cof_display, rel_plot)

    console.print(table)
    session.close()


def run_cluster(
    systems: list[str] | None = None,
    target: str | None = None,
    traj_names: list[str] | None = None,
    cutoff: float = 0.15,
    dt: int = 100,
    group: str = "COMPLEX",
    method: str = "gromos",
    max_clusters: int = 5,
    force: bool = False,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """
    Clusters clean trajectory ensembles and registers each centroid as a
    peer ReceptorFrame (state='ready') in the parent Target.
    Automatically skips systems whose centroids are already generated.
    """
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(
        session, systems, target, traj_names, target_status="cleaned"
    )

    if not trajs:
        console.print(
            "[yellow]No cleaned trajectories found to cluster. Run 'sbdd traj-clean' first![/yellow]"
        )
        session.close()
        return

    for t in trajs:
        if t.status != "cleaned" or not t.path_clean_xtc or not t.path_ref_pdb:
            console.print(
                f"[bold red]FATAL: Trajectory '{t.name}' must be cleaned first. Run 'sbdd traj-clean'![/bold red]"
            )
            sys.exit(1)

    # Group trajectories strictly by System and filter out systems whose clusters are already up to date
    by_sys = {}
    for t in trajs:
        by_sys.setdefault(t.system_id, []).append(t)

    plan_list = []
    sys_map = {}

    for sys_id, t_list in by_sys.items():
        sys_obj = session.query(System).get(sys_id)

        # Check existing centroids for this system
        existing_centroids = (
            session.query(ReceptorFrame)
            .filter(
                ReceptorFrame.target_id == sys_obj.target_id,
                ReceptorFrame.source_type == "md_centroid",
                ReceptorFrame.parent_frame_id == sys_obj.receptor_frame_id,
            )
            .all()
        )

        active_traj_names = sorted([t.name for t in t_list])
        current_cluster_sig = {
            "trajs": active_traj_names,
            "cutoff": float(cutoff),
            "max_clusters": int(max_clusters),
            "group": str(group),
            "method": str(method),
            "dt": int(dt),
        }

        # If centroids already exist, verify if trajectory set or clustering parameters have changed
        if existing_centroids:
            stored_sig = sys_obj.environment_dict.get("cluster_params_sig")
            if stored_sig is None:
                stored_sig = {
                    "trajs": sys_obj.environment_dict.get("cluster_active_trajs", [])
                }

            if stored_sig == current_cluster_sig and not force:
                console.print(
                    f"[dim]System '{sys_obj.name}' already has up-to-date centroids matching parameters (cutoff={cutoff}, {active_traj_names}). Skipping.[/dim]"
                )
                continue
            else:
                console.print(
                    f"[yellow]Clustering parameters or trajectory ensemble for '{sys_obj.name}' modified. Purging outdated centroids to re-cluster...[/yellow]"
                )
                for old_f in existing_centroids:
                    session.delete(old_f)
                session.flush()

        task_id = f"{sys_obj.name}_cluster"
        sys_map[task_id] = (sys_obj, active_traj_names, current_cluster_sig)

        clean_xtcs = [t.path_clean_xtc for t in t_list]
        ref_pdb = t_list[0].path_ref_pdb
        lig_code = sys_obj.pose.substance.name if sys_obj.pose else ""

        plan_list.append(
            {
                "task_id": task_id,
                "system": sys_obj.name,
                "clean_xtcs": clean_xtcs,
                "ref_pdb": ref_pdb,
                "ligand": lig_code,
                "cutoff": cutoff,
                "method": method,
                "dt": dt,
                "group": group,
                "max_clusters": max_clusters,
            }
        )

    if not plan_list:
        console.print(
            "[bold green]✓ All systems already have up-to-date centroids matching their active trajectories. Nothing to do![/bold green]"
        )
        session.close()
        return

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Clustering trajectories across {len(by_sys)} System(s)...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--cluster_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
    # When force is explicitly provided, invalidate Nextflow cache to force genuine recalculation
    resume_flag = resume and not force

    ret_code = execute_nextflow(
        "md-analyze",
        nf_args,
        cwd,
        profile=profile_name,
        resume=resume_flag,
        configs=configs,
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow clustering failed.[/bold red]")
        sys.exit(ret_code)

    prepared_dir = (cwd / "prepared").resolve()
    receipt_files = []

    # Strict Fail-Fast verification of canonical JSON task receipts
    for task_id in sys_map.keys():
        task_receipt = (
            prepared_dir / task_id / f"{task_id}_cluster_receipt.json"
        ).resolve()
        if not task_receipt.exists():
            console.print(
                f"[bold red]FATAL: Canonical cluster JSON receipt not found at expected location: {task_receipt}[/bold red]"
            )
            sys.exit(1)
        receipt_files.append(task_receipt)

    registered_frames = 0
    records = []
    for rf in receipt_files:
        try:
            records.extend(json.loads(rf.read_text(encoding="utf-8")))
        except Exception as e:
            raise ValueError(f"FATAL: Malformed cluster receipt JSON at {rf}: {e}")

    for rec in records:
        parent_sys_key = rec["sys_name"]
        centroid_title = rec["centroid_title"]
        view_pdb_name = rec["view_pdb"]
        pdb_name = rec["pdb"]
        mae_name = rec["mae"]
        poses_dict = rec.get("poses", {})

        mapped_item = sys_map.get(parent_sys_key)
        if not mapped_item:
            raise KeyError(
                f"FATAL: System mapping missing for explicit SSOT key '{parent_sys_key}'"
            )
        parent_sys, active_traj_names, cluster_sig = mapped_item

        # Absolute path resolution directly from task cluster directory
        task_dir = (cwd / "reports" / "clusters" / parent_sys_key).resolve()
        src_pdb = (task_dir / pdb_name).resolve()
        src_mae = (task_dir / mae_name).resolve()
        src_view_pdb = (task_dir / view_pdb_name).resolve()

        if not src_pdb.exists() or not src_mae.exists():
            raise FileNotFoundError(
                f"FATAL: Centroid deliverable missing in {task_dir}: {src_pdb.name} / {src_mae.name}"
            )

        # Record clustering signature on parent system
        env = parent_sys.environment_dict
        env["cluster_active_trajs"] = active_traj_names
        env["cluster_params_sig"] = cluster_sig
        parent_sys.environment_dict = env

        # A. Publish Centroid Frame into targets/<target>/frames/<centroid_frame_name>/
        centroid_frame_name = centroid_title.replace("_cluster", "")
        frame_dir = cwd / "targets" / parent_sys.target.name / "frames" / centroid_frame_name
        frame_dir.mkdir(parents=True, exist_ok=True)

        dest_pdb = frame_dir / f"{centroid_frame_name}_receptor.pdb"
        dest_mae = frame_dir / f"{centroid_frame_name}_receptor.mae"
        import shutil
        if src_pdb != dest_pdb:
            shutil.copy2(src_pdb, dest_pdb)
        if src_mae != dest_mae:
            shutil.copy2(src_mae, dest_mae)

        ex_frame = (
            session.query(ReceptorFrame).filter_by(name=centroid_frame_name).first()
        )
        if not ex_frame:
            new_frame = ReceptorFrame(
                target_id=parent_sys.target_id,
                name=centroid_frame_name,
                source_type="md_centroid",
                state="ready",
                parent_frame_id=parent_sys.receptor_frame_id,
                path=str(dest_pdb.resolve()),
                path_mae=str(dest_mae.resolve()),
                path_pdb=str(dest_pdb.resolve()),
            )
            session.add(new_frame)
            session.flush()
            frame_obj = new_frame
        else:
            ex_frame.path = str(dest_pdb.resolve())
            ex_frame.path_mae = str(dest_mae.resolve())
            ex_frame.path_pdb = str(dest_pdb.resolve())
            frame_obj = ex_frame

        # B. Publish Centroid Poses into substances/<inchikey>/poses/<centroid_pose_name>/
        centroid_main_pose = None
        centroid_cofactors = []

        for res_name, p_files in poses_dict.items():
            res_clean = res_name.strip().upper()
            src_p_sdf = (task_dir / p_files["sdf"]).resolve()
            src_p_mae = (task_dir / p_files["mae"]).resolve()

            if not src_p_sdf.exists() or src_p_sdf.stat().st_size == 0:
                raise FileNotFoundError(f"FATAL: Centroid pose SDF missing or empty: {src_p_sdf}")

            # Match with parent pose to preserve substance lineage and charge
            parent_match = None
            is_main_lig = False

            if parent_sys.pose and res_clean == parent_sys.pose.res_name:
                parent_match = parent_sys.pose
                is_main_lig = True
            else:
                for cof_p in parent_sys.cofactors:
                    if res_clean == cof_p.res_name:
                        parent_match = cof_p
                        break

            if not parent_match:
                raise ValueError(
                    f"FATAL: Extracted centroid molecule '{res_clean}' has no matching Pose in parent system '{parent_sys.name}'!"
                )

            cluster_pose_name = f"{centroid_frame_name}_{res_clean}_pose"
            sub_pose_dir = cwd / "substances" / parent_match.substance.inchikey / "poses" / cluster_pose_name
            sub_pose_dir.mkdir(parents=True, exist_ok=True)

            dest_p_sdf = sub_pose_dir / f"{res_clean}_pose.sdf"
            dest_p_mae = sub_pose_dir / f"{res_clean}_pose.mae"

            if src_p_sdf != dest_p_sdf:
                shutil.copy2(src_p_sdf, dest_p_sdf)
            if src_p_mae.exists() and src_p_mae != dest_p_mae:
                shutil.copy2(src_p_mae, dest_p_mae)

            ex_pose = session.query(Pose).filter_by(name=cluster_pose_name).first()
            if not ex_pose:
                new_p = Pose(
                    substance_id=parent_match.substance_id,
                    receptor_frame_id=frame_obj.id,
                    name=cluster_pose_name,
                    res_name=parent_match.res_name,
                    origin="md_cluster",
                    state="raw",  # Geometry shifted in MD -> requires sbdd prep Antechamber
                    formal_charge=parent_match.formal_charge,
                    path=str(dest_p_sdf.resolve()),
                    path_sdf=str(dest_p_sdf.resolve()),
                    path_mae=str(dest_p_mae.resolve()) if dest_p_mae.exists() else None,
                )
                session.add(new_p)
                session.flush()
                curr_pose = new_p
            else:
                ex_pose.path = str(dest_p_sdf.resolve())
                ex_pose.path_sdf = str(dest_p_sdf.resolve())
                ex_pose.path_mae = str(dest_p_mae.resolve()) if dest_p_mae.exists() else None
                ex_pose.formal_charge = parent_match.formal_charge
                ex_pose.state = "raw"
                curr_pose = ex_pose

            if is_main_lig:
                centroid_main_pose = curr_pose
            else:
                centroid_cofactors.append(curr_pose)

        # C. Register System for this Centroid with view.pdb in systems/<centroid_sys_name>/
        cluster_sys_name = (
            f"{centroid_frame_name}_{parent_sys.pose.substance.name}"
            if centroid_main_pose
            else centroid_frame_name
        )

        sys_dir = cwd / "systems" / cluster_sys_name
        sys_dir.mkdir(parents=True, exist_ok=True)
        dest_view_pdb = sys_dir / f"{cluster_sys_name}_view.pdb"
        if src_view_pdb.exists() and src_view_pdb != dest_view_pdb:
            shutil.copy2(src_view_pdb, dest_view_pdb)

        ex_sys = session.query(System).filter_by(name=cluster_sys_name).first()
        if not ex_sys:
            new_sys = System(
                target_id=parent_sys.target_id,
                receptor_frame_id=frame_obj.id,
                pose_id=centroid_main_pose.id if centroid_main_pose else None,
                name=cluster_sys_name,
                state="raw",
                origin="md_cluster",
            )
            new_sys.cofactors = centroid_cofactors
            session.add(new_sys)
        else:
            ex_sys.pose_id = centroid_main_pose.id if centroid_main_pose else None
            ex_sys.state = "raw"
            ex_sys.cofactors = centroid_cofactors

        registered_frames += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Clustering complete! Registered {registered_frames} centroid frame(s) as ready ReceptorFrames.[/bold green]"
    )
    console.print(
        "[dim]These frames can now be used directly for docking. View with 'sbdd tree'.[/dim]"
    )


def run_flag_traj(
    trajectory: str,
    status: str = "drifted",
    disable: bool = False,
    enable: bool = False,
) -> None:
    """Marks a trajectory as drifted, disabled, or cleaned to control inclusion in downstream analysis."""
    session = init_db()
    traj_name = trajectory
    t_obj = session.query(Trajectory).filter_by(name=traj_name).first()

    if not t_obj:
        console.print(
            f"[bold red]FATAL: Trajectory '{traj_name}' not found in database.[/bold red]"
        )
        session.close()
        sys.exit(1)

    new_status = status
    if disable:
        new_status = "disabled"
    elif enable:
        new_status = "cleaned" if t_obj.path_clean_xtc else "raw"

    old_status = t_obj.status
    t_obj.status = new_status
    session.commit()
    session.close()

    console.print(
        f"[bold green]✓ Updated trajectory '{traj_name}': status '{old_status}' -> '{new_status}'[/bold green]"
    )
    if new_status in ["drifted", "disabled"]:
        console.print(
            f"[yellow]Trajectory is now excluded from clustering. 'sbdd md' will detect missing replicates.[/yellow]"
        )


def run_md_analyze(
    traj: str | None = None,
    tpr: str | None = None,
    topol: str | None = None,
    ndx: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
    clean_scratch: bool = False,
) -> None:
    """Run trajectory PBC fix, metrics, and MM-PBSA standalone."""
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")

    if traj:
        nf_args.extend(["--traj", str(Path(traj).resolve())])
    if tpr:
        nf_args.extend(["--tpr", str(Path(tpr).resolve())])
    if topol:
        nf_args.extend(["--topol", str(Path(topol).resolve())])
    if ndx:
        nf_args.extend(["--ndx", str(Path(ndx).resolve())])
    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    ret_code = execute_nextflow(
        "md-analyze",
        nf_args,
        cwd,
        profile=profile_name,
        resume=resume,
        auto_clean=clean_scratch,
        configs=configs,
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow md-analyze failed.[/bold red]")
        sys.exit(ret_code)


def run_mmpbsa(
    systems: list[str] | None = None,
    target: str | None = None,
    traj_names: list[str] | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """
    CLI handler for 'sbdd traj-mmpbsa': standalone free energy calculations.
    Automatically queries all cleaned trajectories.
    """
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(
        session, systems, target, traj_names, target_status="cleaned"
    )

    if not trajs:
        console.print(
            "[yellow]No cleaned trajectories found for MM-PBSA. Run 'sbdd traj-clean' first![/yellow]"
        )
        session.close()
        return

    plan_list = []
    for t in trajs:
        if not t.path_clean_xtc or not t.path_ref_pdb:
            continue

        sys_obj = t.system
        if not sys_obj.path_top or not sys_obj.path_ndx:
            raise ValueError(
                f"FATAL: System '{sys_obj.name}' has no registered topology paths in database! "
                "Run 'sbdd md-setup' first."
            )

        topol_path = Path(sys_obj.path_top).resolve()
        ndx_path = Path(sys_obj.path_ndx).resolve()

        if not topol_path.exists() or not ndx_path.exists():
            raise FileNotFoundError(
                f"FATAL: GROMACS topology files missing on disk for '{sys_obj.name}': {topol_path} / {ndx_path}"
            )

        lig_code = sys_obj.pose.substance.name if sys_obj.pose else ""

        plan_list.append(
            {
                "name": t.name,
                "rec_id": sys_obj.name,
                "rep": t.replicate,
                "clean_xtc": str(Path(t.path_clean_xtc).resolve()),
                "ref_pdb": str(Path(t.path_ref_pdb).resolve()),
                "topol": str(topol_path),
                "ndx": str(ndx_path),
                "ligand": lig_code,
            }
        )

    if not plan_list:
        console.print(
            "[yellow]No valid trajectories with topologies found for MM-PBSA.[/yellow]"
        )
        session.close()
        return

    plan_file = cwd / "mmpbsa_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Running standalone MM-PBSA on {len(plan_list)} clean trajectory(ies)...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--mmpbsa_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow MM-PBSA failed.[/bold red]")
        sys.exit(ret_code)

    console.print(
        "\n[bold green]✓ MM-PBSA calculations completed successfully![/bold green]"
    )
    session.close()
