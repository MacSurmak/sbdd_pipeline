"""
File: cli/commands/traj.py
Description: Clean 3NF CLI handlers for trajectory solvent stripping, metrics, and clustering.
Enforces stream cleaning in worker scratch without raw file duplication on disk.
"""

import json
import sys
import shutil
import re
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from rich.console import Console
from rich.table import Table
from cli.core.db import (
    get_session,
    init_db,
    Target,
    ReceptorFrame,
    Substance,
    SubstanceState,
    Pose,
    System,
    Trajectory,
    ClusteringRun,
    TrajectoryCluster,
    MmpbsaRun,
)
from cli.core.runner import execute_nextflow
from cli.core.config_loader import (
    load_project_configs,
    build_nextflow_args,
    validate_clustering_config,
    validate_analysis_sampling_config,
)

console = Console()

FORBIDDEN_SCRATCH_TOKENS = [
    "/work/",
    ".nextflow",
    "/tmp/",
    "slurm_tmp",
    "nfs_scratch",
    "nxf.",
]


def assert_permanent_path(path: Path | str, context: str) -> Path:
    """
    Absolute assertion barrier preventing ephemeral Nextflow or Slurm scratch paths
    from entering the CLI runtime, manifests, or SQLite database.
    """
    p = Path(path).resolve()
    p_str = str(p)
    for token in FORBIDDEN_SCRATCH_TOKENS:
        if token in p_str:
            raise AssertionError(
                f"FATAL SCRATCH LEAKAGE DETECTED in {context}: '{p_str}' contains forbidden ephemeral token '{token}'! "
                "Only validated, permanent deliverables located in /shared may be processed."
            )
    return p


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
    session = get_session()
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
        lig_code = t.system.pose.res_name if (t.system and t.system.pose) else ""
        plan_list.append(
            {
                "name": t.name,
                "rec_id": t.system.name,
                "rep": t.replicate,
                "xtcs": [str(c.resolve()) for c in chunks],
                "tpr": str(Path(t.path_tpr).resolve()),
                "ligand": lig_code,
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
    session = get_session()
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
        lig_code = (
            sys_obj.pose.substance.name
            if (sys_obj.pose and sys_obj.pose.substance)
            else ""
        )
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


def save_clustering_to_db(
    session,
    sys_obj: System,
    manifest_path: Path,
    custom_name: str | None = None,
    overwrite: bool = False,
) -> ClusteringRun:
    """
    Atomically registers pure analytical clustering run and its macro-basins into SQLite 3NF.
    Strictly preserves semantic names per system without prepending system name prefixes.
    """
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    task_dir = manifest_path.parent
    eff_align = manifest["align_group"]
    eff_cluster = manifest["cluster_group"]
    eff_dim = manifest["dim_reduction"]
    if "algorithm" not in manifest:
        raise KeyError(
            f"FATAL: Required 'algorithm' key missing in manifest {manifest_path}! "
            "Fallbacks are strictly prohibited by Fail-Fast law."
        )
    eff_algo = manifest["algorithm"]

    configs = load_project_configs(Path.cwd())
    val_sampling = validate_analysis_sampling_config(configs)
    start_ns = float(val_sampling["start_ns"])
    interval_ps = float(val_sampling["sample_interval_ps"])

    clean_group_tag = re.sub(r"[^a-zA-Z0-9_]", "_", eff_cluster).strip("_")
    run_name = custom_name.strip() if custom_name else f"{clean_group_tag}_{eff_algo}"

    labels_p = assert_permanent_path(
        task_dir / "cluster_labels.npy", "Cluster Labels Array"
    )
    time_map_p = assert_permanent_path(
        task_dir / "frame_time_map.json", "Frame Time Map JSON"
    )

    if not labels_p.exists() or not time_map_p.exists():
        session.close()
        raise FileNotFoundError(
            f"FATAL: Required clustering deliverables missing in {task_dir} (labels or time map)!"
        )

    ex_run = session.query(ClusteringRun).filter_by(system_id=sys_obj.id, name=run_name).first()
    if ex_run:
        if not overwrite:
            session.close()
            raise ValueError(
                f"FATAL: ClusteringRun '{run_name}' already exists for system '{sys_obj.name}'! "
                "Use --overwrite to replace it."
            )
        session.delete(ex_run)
        session.flush()

    c_run = ClusteringRun(
        system_id=sys_obj.id,
        name=run_name,
        feature_group=eff_cluster,
        align_group=eff_align,
        algorithm=eff_algo,
        dim_reduction=eff_dim,
        start_ns=start_ns,
        sample_interval_ps=interval_ps,
        parameters_dict=manifest.get("feature_metadata", {}),
        path_labels=str(labels_p.resolve()),
        path_time_map=str(time_map_p.resolve()),
    )
    session.add(c_run)
    session.flush()

    basins = manifest.get("basins", [])
    for b in basins:
        tc = TrajectoryCluster(
            clustering_run_id=c_run.id,
            cluster_index=int(b["raw_cluster_id"]),
            cluster_rank=int(b["rank"]),
            population_fraction=round(float(b["population_pct"]) / 100.0, 4),
            frame_count=int(b["population_count"]),
            medoid_replicate=str(b["medoid_replicate"]),
            medoid_frame=int(b["medoid_frame"]),
            medoid_time_ns=float(b["medoid_time_ns"]),
            compactness_spread=float(b["compactness_spread"]),
        )
        session.add(tc)

    session.commit()
    console.print(
        f"[bold green]✓ Saved analytical ClusteringRun '{run_name}' with {len(basins)} basin(s) into SQLite 3NF.[/bold green]"
    )
    return c_run


def run_extract_medoids(
    clustering: str,
    system: str | None = None,
    ranks: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """
    Dedicated on-demand extractor: orchestrates medoid complex extraction on Slurm cluster nodes.
    Follows Plan & Receipt pattern: CLI constructs plan -> Nextflow runs in $SLURM_TMPDIR ->
    Receipt synchronizes 3NF SQLite database. Master node does ZERO heavy I/O.
    """
    cwd = Path.cwd()
    session = get_session()

    query = session.query(ClusteringRun).filter_by(name=clustering.strip())
    if system:
        sys_obj = session.query(System).filter_by(name=system.strip()).first()
        if not sys_obj:
            session.close()
            raise ValueError(f"FATAL: System '{system}' not found in database!")
        query = query.filter_by(system_id=sys_obj.id)

    matched_runs = query.all()
    if not matched_runs:
        avail_runs = sorted(list(set(r.name for r in session.query(ClusteringRun.name).all())))
        session.close()
        raise ValueError(
            f"FATAL: ClusteringRun '{clustering}' not found in database!\n"
            f"Available clustering runs in project: {avail_runs}"
        )

    plan_tasks = []

    for c_run in matched_runs:
        sys_obj = c_run.system
        target_obj = sys_obj.target
        parent_frame = sys_obj.receptor_frame
        parent_pose = sys_obj.pose

        if ranks:
            try:
                target_ranks = [
                    int(x.strip().replace("#", "")) for x in ranks.split(",") if x.strip()
                ]
            except ValueError:
                session.close()
                raise ValueError(
                    f"FATAL: Invalid ranks specified: '{ranks}'. Expected comma-separated integers."
                )
        else:
            target_ranks = [tc.cluster_rank for tc in c_run.clusters]

        selected_clusters = [tc for tc in c_run.clusters if tc.cluster_rank in target_ranks]
        if not selected_clusters:
            continue

        all_trajs = (
            session.query(Trajectory)
            .filter_by(system_id=sys_obj.id, status="cleaned")
            .all()
        )
        if not all_trajs:
            session.close()
            raise FileNotFoundError(
                f"FATAL: No cleaned trajectories found for system '{sys_obj.name}'!"
            )

        ref_pdb_path = Path(all_trajs[0].path_ref_pdb).resolve()
        clean_xtcs = [str(Path(t.path_clean_xtc).resolve()) for t in all_trajs]

        clusters_spec = [
            {
                "rank": tc.cluster_rank,
                "medoid_frame": tc.medoid_frame,
                "medoid_replicate": tc.medoid_replicate,
                "medoid_time_ns": tc.medoid_time_ns,
                "population": tc.population_fraction,
                "spread": tc.compactness_spread,
            }
            for tc in selected_clusters
        ]

        lig_spec = None
        if parent_pose:
            lig_spec = {
                "res_name": parent_pose.res_name,
                "inchikey": parent_pose.substance.inchikey,
                "template_sdf": str(Path(parent_pose.substance_state.path_sdf).resolve()),
                "state_id": parent_pose.substance_state_id,
            }

        cofactors_spec = [
            {
                "res_name": c.res_name,
                "inchikey": c.substance.inchikey,
                "template_sdf": str(Path(c.substance_state.path_sdf).resolve()),
                "state_id": c.substance_state_id,
                "is_ion": (c.substance.substance_type == "ion"),
            }
            for c in sys_obj.cofactors
        ]

        plan_tasks.append(
            {
                "task_id": f"{sys_obj.name}_{c_run.name}",
                "system_name": sys_obj.name,
                "target_name": target_obj.name,
                "clustering_name": c_run.name,
                "parent_frame_name": parent_frame.name,
                "parent_frame_id": parent_frame.id,
                "feature_group": c_run.feature_group,
                "align_group": c_run.align_group,
                "dim_reduction": c_run.dim_reduction,
                "ref_pdb": str(ref_pdb_path),
                "clean_xtcs": clean_xtcs,
                "clusters": clusters_spec,
                "ligand": lig_spec,
                "cofactors": cofactors_spec,
            }
        )

    if not plan_tasks:
        session.close()
        console.print("[yellow]No matching systems or clusters found to extract.[/yellow]")
        return

    plan_file = cwd / "extract_medoids_plan.json"
    plan_file.write_text(json.dumps(plan_tasks, indent=2), encoding="utf-8")

    receipt_file = cwd / "extract_medoids_receipt.tsv"
    receipt_file.unlink(missing_ok=True)

    console.print(
        f"[bold cyan]Dispatching extraction of medoid complexes across {len(plan_tasks)} system(s) to Slurm cluster...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--extract_medoids_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        session.close()
        console.print("[bold red]Nextflow extraction failed on cluster.[/bold red]")
        sys.exit(ret_code)

    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not emitted by Nextflow. "
            "Extraction failed without generating receipt."
        )

    # Atomic 3NF Database Synchronization based on verified receipt
    frames_created = 0
    poses_created = 0
    systems_created = 0

    frame_lines = []
    pose_lines = []
    system_lines = []

    with open(receipt_file, "r", encoding="utf-8") as rf:
        for line in rf:
            line_str = line.strip()
            if not line_str:
                continue
            parts = line_str.split("\t")
            record_type = parts[0]
            if record_type == "FRAME":
                frame_lines.append(parts)
            elif record_type == "POSE":
                pose_lines.append(parts)
            elif record_type == "SYSTEM":
                system_lines.append(parts)

    # Phase 1: Register all ReceptorFrames and flush to populate primary integer keys
    frame_cache = {}
    for parts in frame_lines:
        _, t_name, f_name, rel_pdb, rel_mae, p_id_str, r_str, pop_str, sp_str, f_grp, a_grp, d_red = parts
        abs_pdb = assert_permanent_path(cwd / rel_pdb, f"Receptor PDB {f_name}")
        abs_mae = assert_permanent_path(cwd / rel_mae, f"Receptor MAE {f_name}") if rel_mae != "None" else None

        t_obj = session.query(Target).filter_by(name=t_name).first()
        ex_f = session.query(ReceptorFrame).filter_by(name=f_name).first()
        if not ex_f:
            new_frame = ReceptorFrame(
                target_id=t_obj.id,
                name=f_name,
                source_type="md_centroid",
                parent_frame_id=int(p_id_str) if p_id_str != "None" else None,
                feature_group=f_grp,
                align_group=a_grp,
                dim_reduction_method=d_red,
                cluster_rank=int(r_str),
                population=float(pop_str),
                compactness_spread=float(sp_str),
                state="ready",
                path=str(abs_pdb),
                path_pdb=str(abs_pdb),
                path_mae=str(abs_mae) if abs_mae else None,
            )
            session.add(new_frame)
            frame_cache[f_name] = new_frame
            frames_created += 1
        else:
            ex_f.path = str(abs_pdb)
            ex_f.path_pdb = str(abs_pdb)
            if abs_mae:
                ex_f.path_mae = str(abs_mae)
            ex_f.state = "ready"
            ex_f.population = float(pop_str)
            ex_f.compactness_spread = float(sp_str)
            frame_cache[f_name] = ex_f

    session.flush()

    # Phase 2: Register all Poses linked to guaranteed ReceptorFrame IDs
    pose_cache = {}
    for parts in pose_lines:
        _, f_name, p_name, inchikey, res_name, st_id_str, rel_sdf, rel_mae, lig_flag = parts
        abs_sdf = assert_permanent_path(cwd / rel_sdf, f"Pose SDF {p_name}")
        abs_mae = assert_permanent_path(cwd / rel_mae, f"Pose MAE {p_name}") if rel_mae != "None" else None

        frame_obj = frame_cache.get(f_name) or session.query(ReceptorFrame).filter_by(name=f_name).first()
        if not frame_obj:
            session.close()
            raise ValueError(f"FATAL: ReceptorFrame '{f_name}' referenced by pose '{p_name}' not found!")

        state_obj = session.query(SubstanceState).get(int(st_id_str))
        if not state_obj:
            session.close()
            raise ValueError(f"FATAL: SubstanceState ID {st_id_str} referenced by pose '{p_name}' not found!")

        ex_p = session.query(Pose).filter_by(name=p_name).first()
        if not ex_p:
            new_p = Pose(
                substance_state_id=state_obj.id,
                receptor_frame_id=frame_obj.id,
                name=p_name,
                res_name=res_name,
                origin="md_centroid",
                state="ready",
                path=str(abs_sdf),
                path_sdf=str(abs_sdf),
                path_mae=str(abs_mae) if abs_mae else None,
            )
            session.add(new_p)
            pose_cache[p_name] = new_p
            poses_created += 1
        else:
            ex_p.path = str(abs_sdf)
            ex_p.path_sdf = str(abs_sdf)
            if abs_mae:
                ex_p.path_mae = str(abs_mae)
            ex_p.state = "ready"
            pose_cache[p_name] = ex_p

    session.flush()

    # Phase 3: Assemble physical Systems binding Frames, Ligand Pose, and Cofactor Poses
    for parts in system_lines:
        _, t_name, s_name, f_name, lig_pname, cof_names_str = parts
        t_obj = session.query(Target).filter_by(name=t_name).first()
        frame_obj = frame_cache.get(f_name) or session.query(ReceptorFrame).filter_by(name=f_name).first()

        lig_pose = pose_cache.get(lig_pname) or (session.query(Pose).filter_by(name=lig_pname).first() if lig_pname else None)
        cofactor_poses = []
        if cof_names_str:
            for cp in cof_names_str.split(","):
                c_name_clean = cp.strip()
                if c_name_clean:
                    found_cof = pose_cache.get(c_name_clean) or session.query(Pose).filter_by(name=c_name_clean).first()
                    if found_cof:
                        cofactor_poses.append(found_cof)

        ex_sys = session.query(System).filter_by(name=s_name).first()
        if not ex_sys:
            new_sys = System(
                target_id=t_obj.id,
                receptor_frame_id=frame_obj.id,
                pose_id=lig_pose.id if lig_pose else None,
                name=s_name,
                state="raw",
                origin="md_centroid",
            )
            new_sys.cofactors = cofactor_poses
            session.add(new_sys)
            systems_created += 1
        else:
            ex_sys.pose_id = lig_pose.id if lig_pose else None
            ex_sys.cofactors = cofactor_poses
            ex_sys.state = "raw"

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Cluster extraction completed via Slurm! Synchronized with SQLite 3NF:[/bold green]\n"
        f"  - Receptor Frames: {frames_created} created/updated (PDB v3.3 + MAE)\n"
        f"  - Small Molecule Poses: {poses_created} created/updated\n"
        f"  - Complete Physical Systems: {systems_created} created/updated\n"
    )


def run_cluster(
    systems: list[str] | None = None,
    target: str | None = None,
    traj_names: list[str] | None = None,
    align_group: str | None = None,
    cluster_group: str | None = None,
    dim_reduction: str | None = None,
    algorithm: str | None = None,
    n_clusters: int | None = None,
    rmsd_cutoff: float | None = None,
    min_population: float | None = None,
    name: str | None = None,
    save: bool = False,
    overwrite: bool = False,
    auto: bool = False,
    force: bool = False,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """
    Orchestrates conformational macro-basin clustering via PCA/TICA + UMAP and K-Means medoids.
    Presents an interactive Rich UI summary table for visual validation before promoting
    centroids into SQLite 3NF as ready ReceptorFrame entities (Lineage Variant 1).
    """
    cwd = Path.cwd()
    session = get_session()
    configs = load_project_configs(cwd)
    cluster_cfg = validate_clustering_config(configs)

    # Apply explicit CLI overrides onto configuration dictionary before Nextflow argument construction
    eff_align = align_group if align_group is not None else cluster_cfg["align_group"]
    eff_cluster = (
        cluster_group if cluster_group is not None else cluster_cfg["cluster_group"]
    )
    eff_dim = (
        dim_reduction if dim_reduction is not None else cluster_cfg["dim_reduction"]
    )
    eff_algo = (
        algorithm if algorithm is not None else cluster_cfg.get("algorithm", "hdbscan")
    )
    eff_k = n_clusters if n_clusters is not None else cluster_cfg.get("n_clusters", 3)
    eff_cutoff = (
        rmsd_cutoff if rmsd_cutoff is not None else cluster_cfg.get("rmsd_cutoff", 1.5)
    )
    eff_min_pop = (
        min_population
        if min_population is not None
        else cluster_cfg["min_cluster_population"]
    )

    configs.setdefault("md_simulation", {}).setdefault("clustering", {})[
        "align_group"
    ] = eff_align
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})[
        "cluster_group"
    ] = eff_cluster
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})[
        "dim_reduction"
    ] = eff_dim
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})[
        "algorithm"
    ] = eff_algo
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})[
        "rmsd_cutoff"
    ] = eff_cutoff
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})[
        "n_clusters"
    ] = eff_k
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})[
        "min_cluster_population"
    ] = eff_min_pop

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
            session.close()
            sys.exit(1)

    by_sys = {}
    for t in trajs:
        by_sys.setdefault(t.system_id, []).append(t)

    plan_list = []
    sys_lookup = {}

    for sys_id, t_list in by_sys.items():
        sys_obj = session.query(System).get(sys_id)
        sys_lookup[sys_obj.name] = (sys_obj, t_list)

        clean_xtcs = [str(Path(t.path_clean_xtc).resolve()) for t in t_list]
        rep_names = [f"rep_{t.replicate}" for t in t_list]
        ref_pdb = str(Path(t_list[0].path_ref_pdb).resolve())
        ndx_path = (
            str(Path(sys_obj.path_ndx).resolve())
            if sys_obj.path_ndx and Path(sys_obj.path_ndx).exists()
            else ""
        )
        lig_code = sys_obj.pose.res_name if sys_obj.pose else ""

        if eff_cluster.lower() == "pocket" and not lig_code:
            console.print(
                f"[bold red]FATAL: System '{sys_obj.name}' has no registered ligand pose! "
                "Cannot resolve pocket without reference ligand. Use --cluster-group Backbone or custom residue range.[/bold red]"
            )
            session.close()
            sys.exit(1)

        task_id = f"{sys_obj.name}_{eff_align}_{eff_cluster}"

        plan_list.append(
            {
                "task_id": task_id,
                "system": sys_obj.name,
                "clean_xtcs": clean_xtcs,
                "replicate_names": rep_names,
                "ref_pdb": ref_pdb,
                "ndx": ndx_path,
                "align_group": eff_align,
                "cluster_group": eff_cluster,
                "ligand": lig_code,
            }
        )

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w", encoding="utf-8") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Dispatching macro-basin clustering (Align: {eff_align}, Features: {eff_cluster}, K={eff_k})...[/bold cyan]"
    )

    receipt_file = cwd / "cluster_receipt.tsv"
    receipt_file.unlink(missing_ok=True)

    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(
        [
            "--cluster_plan",
            str(plan_file.resolve()),
        ]
    )

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
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
        session.close()
        console.print("[bold red]Nextflow clustering execution failed.[/bold red]")
        sys.exit(ret_code)

    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not emitted by Nextflow."
        )

    receipt_records = []
    with open(receipt_file, "r", encoding="utf-8") as rf:
        for line in rf:
            l_str = line.strip()
            if not l_str:
                continue
            parts = l_str.split("\t")
            if len(parts) == 4:
                receipt_records.append(
                    {
                        "system": parts[0],
                        "align_group": parts[1],
                        "cluster_group": parts[2],
                        "manifest_path": Path(parts[3]).resolve(),
                    }
                )

    session.close()

    for rec_meta in receipt_records:
        sys_name = rec_meta["system"]
        m_path = rec_meta["manifest_path"]

        if not m_path.exists() or m_path.stat().st_size == 0:
            raise FileNotFoundError(
                f"FATAL: Published cluster manifest missing on disk: {m_path}"
            )

        with open(m_path, "r", encoding="utf-8") as jf:
            manifest = json.load(jf)

        basins = manifest.get("basins", [])
        if not basins:
            console.print(
                f"[yellow]Warning: No basins found in manifest for system '{sys_name}'.[/yellow]"
            )
            continue

        table = Table(
            title=f"Conformational Macro-Basins: {sys_name} [Align: {eff_align}, Features: {eff_cluster}]",
            show_header=True,
        )
        table.add_column("Selection", justify="center")
        table.add_column("Rank", style="bold cyan", justify="center")
        table.add_column("Population (%)", justify="right")
        table.add_column("Frame / Time (ns)", justify="center")
        table.add_column("Source Replicate", justify="center")
        table.add_column("Compactness (Å)", justify="right")
        table.add_column("Recommendation Verdict", justify="left")

        recommended_ranks = []
        for b in basins:
            rank = int(b["rank"])
            pop = float(b["population_pct"])
            is_rec = pop >= float(eff_min_pop)

            if is_rec:
                recommended_ranks.append(rank)
                sel_box = "[bold green][✓][/bold green]"
                verdict = (
                    "[bold green]Recommended (Major Basin)[/bold green]"
                    if rank == 1
                    else "[green]Recommended[/green]"
                )
            else:
                sel_box = "[dim][ ][/dim]"
                verdict = f"[dim]Low population (< {eff_min_pop}%)[/dim]"

            table.add_row(
                sel_box,
                f"#{rank}",
                f"{pop:.1f}%",
                f"{b['medoid_frame']} ({b['medoid_time_ns']:.1f} ns)",
                str(b["medoid_replicate"]),
                f"{b['compactness_spread']:.2f} Å",
                verdict,
            )

        if not recommended_ranks and basins:
            recommended_ranks.append(1)

        html_report_name = Path(manifest.get("html_report", "")).name
        pml_session_name = Path(manifest.get("pymol_session", "")).name

        html_report_path = assert_permanent_path(
            m_path.parent / html_report_name, "Interactive Landscape HTML"
        )
        pml_session_path = assert_permanent_path(
            m_path.parent / pml_session_name, "Portable PyMOL Session"
        )

        console.print(table)
        console.print(
            f"📊 [bold cyan]Interactive Landscape (FES):[/bold cyan] {html_report_path}"
        )
        console.print(
            f"🔍 [bold cyan]Portable PyMOL 3D Session:[/bold cyan] pymol {pml_session_path}\n"
        )

        should_save = save or auto
        chosen_name = name
        if not should_save and sys.stdin.isatty():
            try:
                user_input = (
                    input(f"Save clustering for '{sys_name}' to database? [y/N]: ")
                    .strip()
                    .lower()
                )
                should_save = user_input in ["y", "yes"]
                if should_save and not chosen_name:
                    def_name = f"{re.sub(r'[^a-zA-Z0-9_]', '_', eff_cluster).strip('_')}_{manifest['algorithm']}"
                    name_input = input(f"Enter identifier name for this clustering run [default: {def_name}]: ").strip()
                    chosen_name = name_input if name_input else def_name
            except (EOFError, KeyboardInterrupt):
                console.print(
                    "\n[yellow]Preview closed. Database not modified.[/yellow]"
                )
                return

        if should_save:
            session_save = get_session()
            sys_db_obj = session_save.query(System).filter_by(name=sys_name).first()
            if sys_db_obj:
                save_clustering_to_db(
                    session=session_save,
                    sys_obj=sys_db_obj,
                    manifest_path=m_path,
                    custom_name=chosen_name,
                    overwrite=overwrite,
                )
            session_save.close()
        else:
            console.print(
                f"[dim]Preview mode: reports preserved in {m_path.parent}. "
                "Database not modified. Pass --save to commit this analysis.[/dim]"
            )


def run_flag_traj(
    trajectory: str,
    status: str = "drifted",
    disable: bool = False,
    enable: bool = False,
) -> None:
    """Marks a trajectory as drifted, disabled, or cleaned to control inclusion in downstream analysis."""
    session = get_session()
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


def run_stratify(
    system: str | None = None,
    all_systems: bool = False,
    clustering: str | None = None,
    temp_k: float | None = None,
) -> None:
    """
    CLI handler for 'sbdd traj-stratify'.
    Projects independent MM-PBSA measurements Delta_G(t) onto conformational macro-basins C(t).
    Calculates basin free energy moments, overall Boltzmann aggregate, and updates SQLite 3NF.
    """
    cwd = Path.cwd()
    session = get_session()
    configs = load_project_configs(cwd)
    cluster_cfg = validate_clustering_config(configs)
    plateau_cfg = cluster_cfg["plateau_detection"]

    target_systems = []
    if system and not all_systems:
        sys_obj = session.query(System).filter_by(name=system).first()
        if not sys_obj:
            session.close()
            raise ValueError(f"FATAL: System '{system}' not found in database!")
        target_systems = [sys_obj]
    else:
        # Default behavior or --all: process all systems that have both MmpbsaRun and ClusteringRun
        target_systems = (
            session.query(System)
            .filter(System.mmpbsa_runs.any(), System.clustering_runs.any())
            .order_by(System.name)
            .all()
        )
        if not target_systems:
            session.close()
            raise FileNotFoundError(
                "FATAL: No physical Systems found having both MM-PBSA and Clustering runs! "
                "Execute 'sbdd traj-mmpbsa' and 'sbdd traj-cluster --save' first."
            )

    console.print(
        f"[bold cyan]Stratifying free energy across {len(target_systems)} physical system(s)...[/bold cyan]\n"
    )

    stratify_bin = (
        Path(__file__).resolve().parent.parent.parent / "bin" / "stratify_mmpbsa.py"
    )
    if not stratify_bin.exists():
        session.close()
        raise FileNotFoundError(
            f"FATAL: Stratification worker not found at {stratify_bin}"
        )

    mmpbsa_reports_dir = cwd / "reports" / "mmpbsa"

    for sys_obj in target_systems:
        m_run = (
            session.query(MmpbsaRun)
            .filter_by(system_id=sys_obj.id)
            .order_by(MmpbsaRun.created_at.desc())
            .first()
        )
        if not m_run:
            console.print(f"[yellow]Skipping '{sys_obj.name}': no MM-PBSA runs registered.[/yellow]")
            continue

        c_runs = (
            session.query(ClusteringRun)
            .filter_by(system_id=sys_obj.id)
            .order_by(ClusteringRun.created_at.desc())
            .all()
        )
        if not c_runs:
            console.print(f"[yellow]Skipping '{sys_obj.name}': no ClusteringRun registered.[/yellow]")
            continue

        if clustering and len(target_systems) == 1:
            target_c_run = next((cr for cr in c_runs if cr.name == clustering.strip()), None)
            if not target_c_run:
                session.close()
                raise ValueError(f"FATAL: ClusteringRun '{clustering}' not found for '{sys_obj.name}'!")
        else:
            target_c_run = c_runs[0]
            console.print(f"[dim]System '{sys_obj.name}': using ClusteringRun '{target_c_run.name}'[/dim]")

        # Provenance verification
        if (
            abs(float(m_run.start_ns) - float(target_c_run.start_ns)) > 1e-3
            or abs(float(m_run.sample_interval_ps) - float(target_c_run.sample_interval_ps)) > 1e-3
        ):
            session.close()
            raise ValueError(
                f"FATAL SSOT PROVENANCE CONFLICT for '{sys_obj.name}': MM-PBSA grid ({m_run.start_ns} ns, {m_run.sample_interval_ps} ps) "
                f"diverges from ClusteringRun ({target_c_run.start_ns} ns, {target_c_run.sample_interval_ps} ps)!"
            )

        task_dir = Path(target_c_run.path_labels).parent
        out_json = task_dir / "stratification_summary.json"
        out_plot = task_dir / f"{target_c_run.name}_stratified_rmsd.png"
        out_ts = task_dir / "stratified_timeseries.tsv"
        effective_temp = float(temp_k if temp_k is not None else m_run.temperature)

        cmd = [
            sys.executable,
            str(stratify_bin.resolve()),
            "--system",
            sys_obj.name,
            "--labels",
            str(target_c_run.path_labels),
            "--time_map",
            str(target_c_run.path_time_map),
            "--mmpbsa_dir",
            str(mmpbsa_reports_dir.resolve()),
            "--metrics_dir",
            str((cwd / "md_runs" / sys_obj.name).resolve()),
            "--temp_k",
            str(effective_temp),
            "--paradigm",
            str(plateau_cfg["paradigm"]),
            "--smoothing_tau_ns",
            str(plateau_cfg["smoothing_tau_ns"]),
            "--significance_z_crit",
            str(plateau_cfg["significance_z_crit"]),
            "--dilation_buffer_ns",
            str(plateau_cfg["dilation_buffer_ns"]),
            "--min_dwell_ns",
            str(plateau_cfg["min_dwell_ns"]),
            "--out_json",
            str(out_json.resolve()),
            "--out_plot",
            str(out_plot.resolve()),
            "--out_timeseries",
            str(out_ts.resolve()),
        ]
        if plateau_cfg["enforce_cluster_continuity"]:
            cmd.append("--enforce_cluster_continuity")

        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            session.close()
            raise RuntimeError(
                f"FATAL: Stratification failed for '{sys_obj.name}':\n{res.stderr}\n{res.stdout}"
            )

        with open(out_json, "r", encoding="utf-8") as f:
            strat_data = json.load(f)

        boltz_total = float(strat_data["boltzmann_total_delta_g"])
        m_run.delta_g_boltzmann_total = boltz_total

        cluster_stats_map = {item["cluster_id"]: item for item in strat_data.get("clusters", [])}
        for tc in target_c_run.clusters:
            if tc.cluster_index in cluster_stats_map:
                c_info = cluster_stats_map[tc.cluster_index]
                tc.delta_g_mean = float(c_info["delta_g_mean"])
                tc.delta_g_sem = float(c_info["delta_g_sem"])
                tc.delta_g_std = float(c_info["delta_g_std"])

        session.commit()

        table = Table(
            title=f"Conformation-Conditioned Free Energy Stratification: {target_c_run.name}",
            show_header=True,
        )
        table.add_column("Basin", style="bold cyan", justify="center")
        table.add_column("Plateau Pts", justify="right")
        table.add_column("Total Pts", justify="right")
        table.add_column("Population", justify="right")
        table.add_column("ΔG Mean ± SEM (kcal/mol)", style="bold green", justify="center")
        table.add_column("Std Dev (kcal/mol)", justify="right")

        for c_info in strat_data.get("clusters", []):
            table.add_row(
                f"Basin #{c_info['cluster_id']}",
                str(c_info["plateau_samples"]),
                str(c_info["total_samples"]),
                f"{c_info['population_fraction']*100:.1f}%",
                f"{c_info['delta_g_mean']:.2f} ± {c_info['delta_g_sem']:.2f}",
                f"{c_info['delta_g_std']:.2f}",
            )

        console.print(table)
        console.print(
            f"═══════════════════════════════════════════════════════════════════════════\n"
            f" [bold magenta]{sys_obj.name} (Boltzmann ΔG_total):[/bold magenta] [bold green]{boltz_total:.2f} kcal/mol[/bold green]\n"
            f"═══════════════════════════════════════════════════════════════════════════"
        )
        console.print(f"📈 [dim]Plot:[/dim] {out_plot}")
        console.print(f"📊 [dim]Timeseries:[/dim] {out_ts}\n")

    session.close()
    console.print("[bold green]✓ All requested systems successfully stratified and synchronized with SQLite 3NF.[/bold green]")


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
    session = get_session()
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
        if not sys_obj.pose:
            session.close()
            raise ValueError(
                f"FATAL: System '{sys_obj.name}' is an Apo system (no ligand pose registered). "
                "MM-PBSA binding free energy requires a Holo system with a bound ligand!"
            )

        if not sys_obj.path_dry_prmtop or not Path(sys_obj.path_dry_prmtop).exists():
            session.close()
            raise FileNotFoundError(
                f"FATAL: System '{sys_obj.name}' is missing dry AMBER topology (path_dry_prmtop)! "
                "Run 'sbdd md-setup --force' to compile dry complex topology."
            )

        dry_prmtop_path = Path(sys_obj.path_dry_prmtop).resolve()
        clean_dir = Path(t.path_clean_xtc).parent
        clean_ndx_path = clean_dir / "complex_clean.ndx"

        if not clean_ndx_path.exists():
            session.close()
            raise FileNotFoundError(
                f"FATAL: Clean index file '{clean_ndx_path}' not found for trajectory '{t.name}'. "
                "Trajectory must be cleaned with updated pipeline. Run 'sbdd traj-clean'."
            )

        lig_code = sys_obj.pose.res_name

        plan_list.append(
            {
                "name": t.name,
                "rec_id": sys_obj.name,
                "rep": t.replicate,
                "clean_xtc": str(Path(t.path_clean_xtc).resolve()),
                "ref_pdb": str(Path(t.path_ref_pdb).resolve()),
                "clean_ndx": str(clean_ndx_path.resolve()),
                "dry_prmtop": str(dry_prmtop_path),
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
        session.close()
        console.print("[bold red]Nextflow MM-PBSA failed.[/bold red]")
        sys.exit(ret_code)

    # Synchronize and register clean DELTA free energies into SQLite 3NF with strict Fail-Fast validation
    reports_dir = cwd / "reports" / "mmpbsa"
    if not reports_dir.exists():
        session.close()
        raise FileNotFoundError(
            f"FATAL: Reports directory not found at {reports_dir}. MM-PBSA output delivery failed."
        )

    # Deduplicate systems: an MmpbsaRun represents the aggregate ensemble of a System, not individual replicates
    unique_sys_names = sorted(list(set(item["rec_id"] for item in plan_list)))
    registered_runs = 0

    for sys_name in unique_sys_names:
        sys_obj = session.query(System).filter_by(name=sys_name).first()
        if not sys_obj:
            session.close()
            raise ValueError(
                f"FATAL: Target System '{sys_name}' declared in plan is not registered in SQLite database!"
            )

        sys_reps = sorted(
            [d for d in reports_dir.glob("*") if sys_name in d.name and d.is_dir()]
        )
        if not sys_reps:
            session.close()
            raise FileNotFoundError(
                f"FATAL: No replicate output directories found for system '{sys_name}' in {reports_dir}!"
            )

        all_frames_dg = []
        rep_means = []

        for r_dir in sys_reps:
            csv_f = r_dir / "thermodynamics_mmpbsa_pb.csv"
            if not csv_f.exists() or csv_f.stat().st_size == 0:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: Missing or empty MM-PBSA deliverable for replicate in {r_dir}!"
                )

            df = pd.read_csv(csv_f)
            target_col = next((c for c in df.columns if "TOTAL" in c.upper()), None)
            if not target_col:
                session.close()
                raise KeyError(
                    f"FATAL: Required 'TOTAL' or 'DELTA TOTAL' column missing in {csv_f}! "
                    f"Available columns: {list(df.columns)}"
                )

            numeric_vals = (
                pd.to_numeric(df[target_col], errors="coerce").dropna().values
            )
            if len(numeric_vals) == 0:
                session.close()
                raise ValueError(
                    f"FATAL: Zero valid numeric energy values extracted from column '{target_col}' in {csv_f}!"
                )

            rep_mean = float(np.mean(numeric_vals))
            rep_means.append(rep_mean)
            all_frames_dg.extend(numeric_vals.tolist())

        if not rep_means:
            session.close()
            raise ValueError(
                f"FATAL: Failed to extract any valid replicate energy means for system '{sys_name}'!"
            )

        sys_mean = float(np.mean(all_frames_dg))
        # Evaluate standard error across replicate means (accounting for between-replicate variance)
        if len(rep_means) > 1:
            sys_sem = float(np.std(rep_means, ddof=1) / np.sqrt(len(rep_means)))
        else:
            sys_sem = float(np.std(all_frames_dg, ddof=1) / np.sqrt(len(all_frames_dg)))

        consolidated_csv = reports_dir / f"{sys_name}_consolidated_mmpbsa.csv"
        pd.DataFrame({"Delta_G": all_frames_dg}).to_csv(consolidated_csv, index=False)

        val_sampling = validate_analysis_sampling_config(configs)
        eff_start_ns = float(val_sampling["start_ns"])
        eff_interval_ps = float(val_sampling["sample_interval_ps"])

        ex_run = session.query(MmpbsaRun).filter_by(system_id=sys_obj.id).first()
        if not ex_run:
            m_run = MmpbsaRun(
                system_id=sys_obj.id,
                start_ns=eff_start_ns,
                sample_interval_ps=eff_interval_ps,
                temperature=float(configs.get("mmpbsa", {}).get("temperature", 310.15)),
                salt_conc=float(configs.get("mmpbsa", {}).get("salt_conc", 0.15)),
                path_csv=str(consolidated_csv.resolve()),
                delta_g_overall_mean=round(sys_mean, 3),
                delta_g_overall_sem=round(sys_sem, 3),
            )
            session.add(m_run)
        else:
            ex_run.start_ns = eff_start_ns
            ex_run.sample_interval_ps = eff_interval_ps
            ex_run.temperature = float(
                configs.get("mmpbsa", {}).get("temperature", 310.15)
            )
            ex_run.salt_conc = float(configs.get("mmpbsa", {}).get("salt_conc", 0.15))
            ex_run.delta_g_overall_mean = round(sys_mean, 3)
            ex_run.delta_g_overall_sem = round(sys_sem, 3)
            ex_run.path_csv = str(consolidated_csv.resolve())

        session.flush()
        registered_runs += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ MM-PBSA calculations completed! Registered thermodynamic profiles for {registered_runs} System(s) in SQLite 3NF.[/bold green]"
    )
