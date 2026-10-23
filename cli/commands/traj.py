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
)


import shutil
from cli.core.runner import execute_nextflow
from cli.core.config_loader import (
    load_project_configs,
    build_nextflow_args,
    validate_clustering_config,
)

console = Console()

FORBIDDEN_SCRATCH_TOKENS = ["/work/", ".nextflow", "/tmp/", "slurm_tmp", "nfs_scratch", "nxf."]


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


def run_cluster_promote(manifest_path: str, ranks: list[int]) -> None:
    """
    Promotes specific cluster centroids from an existing manifest into SQLite 3NF
    without re-running Nextflow or MDAnalysis.
    Guarantees transactional integrity via idempotent Upsert, eliminating UNIQUE constraint crashes.
    """
    cwd = Path.cwd()
    session = get_session()
    m_p = Path(manifest_path).resolve()

    if not m_p.exists() or m_p.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(f"FATAL: Clusters manifest not found at {m_p}!")

    with open(m_p, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    basins = manifest.get("basins", [])
    task_dir = m_p.parent

    # Single Source of Truth: system name is resolved directly from manifest metadata (zero filesystem heuristics)
    sys_name = manifest.get("system")
    if not sys_name:
        session.close()
        raise KeyError(f"FATAL: Canonical 'system' metadata missing in clusters manifest {m_p}!")

    sys_obj = session.query(System).filter_by(name=sys_name).first()
    if not sys_obj:
        session.close()
        raise ValueError(f"FATAL: Parent system '{sys_name}' declared in manifest is not registered in SQLite database!")

    eff_align = manifest.get("align_group")
    eff_cluster = manifest.get("cluster_group")
    eff_dim = manifest.get("dim_reduction", "pca")

    if not eff_align or not eff_cluster:
        session.close()
        raise KeyError(f"FATAL: Missing 'align_group' or 'cluster_group' in manifest {m_p}!")

    # Prefer embedded structures list from clusters_manifest.json; fallback to extracted_manifest.json if needed
    extracted_data = manifest.get("structures")
    if not extracted_data:
        extract_manifest_p = task_dir / "extracted_manifest.json"
        if not extract_manifest_p.exists():
            session.close()
            raise FileNotFoundError(f"FATAL: Extraction manifest missing in {task_dir}!")
        with open(extract_manifest_p, "r", encoding="utf-8") as f:
            extracted_data = json.load(f)

    extracted_by_rank = {int(item["rank"]): item for item in extracted_data}
    all_trajs = session.query(Trajectory).filter_by(system_id=sys_obj.id).all()
    promoted = 0
    selected_ranks_set = set(ranks)

    # Generic Lineage Invariant: query all existing centroids belonging to this exact ensemble lineage.
    # Any centroid whose rank is not in the newly selected set is obsolete and must be pruned.
    stale_centroids = (
        session.query(ReceptorFrame)
        .filter(
            ReceptorFrame.target_id == sys_obj.target_id,
            ReceptorFrame.parent_frame_id == sys_obj.receptor_frame_id,
            ReceptorFrame.feature_group == eff_cluster,
            ReceptorFrame.align_group == eff_align,
            ReceptorFrame.dim_reduction_method == eff_dim,
            ReceptorFrame.source_type == "md_centroid",
            ReceptorFrame.cluster_rank.notin_(selected_ranks_set),
        )
        .all()
    )

    for stale in stale_centroids:
        # Fail-Fast dependency guard: prevent accidental destruction of heavy downstream MD trajectories
        for s in stale.systems:
            if s.trajectories:
                session.close()
                raise ValueError(
                    f"FATAL: Cannot prune obsolete centroid frame '{stale.name}' (rank {stale.cluster_rank}): "
                    f"system '{s.name}' has computed MD trajectories in database! "
                    "Manually resolve or retain this rank."
                )

        # Synchronously purge physical disk deliverables for frame, poses, and systems
        if stale.path:
            f_dir = Path(stale.path).parent
            if f_dir.exists() and f_dir.name == stale.name:
                shutil.rmtree(f_dir, ignore_errors=True)

        for p in stale.poses:
            if p.path:
                p_dir = Path(p.path).parent
                if p_dir.exists() and p_dir.name == p.name:
                    shutil.rmtree(p_dir, ignore_errors=True)

        for s in stale.systems:
            s_dir = cwd / "systems" / s.name
            if s_dir.exists():
                shutil.rmtree(s_dir, ignore_errors=True)

        console.print(
            f"[yellow]Pruned obsolete centroid frame '{stale.name}' (rank {stale.cluster_rank} not in {sorted(list(selected_ranks_set))}) and purged disk artifacts.[/yellow]"
        )
        session.delete(stale)

    session.flush()

    for rank in sorted(list(selected_ranks_set)):
        b_meta = next((b for b in basins if int(b["rank"]) == rank), None)
        ext_meta = extracted_by_rank.get(rank)
        if not b_meta or not ext_meta:
            console.print(f"[yellow]Warning: Basin #{rank} not found in manifest deliverables. Skipping.[/yellow]")
            continue

        raw_rec_name = Path(ext_meta["receptor_pdb"]).name
        raw_view_name = Path(ext_meta["complex_pdb"]).name

        src_rec_pdb = assert_permanent_path(task_dir / raw_rec_name, f"Basin #{rank} Source Receptor PDB")
        src_rec_mae = assert_permanent_path(task_dir / raw_rec_name.replace(".pdb", ".mae"), f"Basin #{rank} Source Receptor MAE")
        src_view_pdb = assert_permanent_path(task_dir / raw_view_name, f"Basin #{rank} Source Complex PDB")

        if not src_rec_pdb.exists() or src_rec_pdb.stat().st_size == 0:
            raise FileNotFoundError(f"FATAL: Centroid receptor PDB missing on disk: {src_rec_pdb}")
        if not src_rec_mae.exists() or src_rec_mae.stat().st_size == 0:
            raise FileNotFoundError(f"FATAL: Centroid receptor MAE missing on disk: {src_rec_mae}")

        # Strict Fail-Fast trajectory replicate matching (zero silent fallbacks to all_trajs[0])
        src_rep_name = str(b_meta["medoid_replicate"])
        matched_traj = next(
            (t for t in all_trajs if t.name == src_rep_name or t.replicate == int(src_rep_name.replace("rep_", ""))),
            None
        )
        if not matched_traj:
            session.close()
            raise ValueError(
                f"FATAL: Source trajectory replicate '{src_rep_name}' for basin #{rank} "
                f"is not registered in database for system '{sys_obj.name}'!"
            )

        # Canonical naming: ReceptorFrame incorporates parent frame name to eliminate UNIQUE collision across multiple models
        parent_frame_name = sys_obj.receptor_frame.name
        new_frame_name = f"{parent_frame_name}_{eff_cluster.lower()}_c{rank}"
        dest_frame_dir = assert_permanent_path(cwd / "targets" / sys_obj.target.name / "frames" / new_frame_name, f"Basin #{rank} Frame Directory")

        # If existing frame directory on disk differs from new canonical directory, purge the obsolete one
        existing_frame_check = session.query(ReceptorFrame).filter(
            ReceptorFrame.target_id == sys_obj.target_id,
            ReceptorFrame.parent_frame_id == sys_obj.receptor_frame_id,
            ReceptorFrame.feature_group == eff_cluster,
            ReceptorFrame.align_group == eff_align,
            ReceptorFrame.dim_reduction_method == eff_dim,
            ReceptorFrame.cluster_rank == rank,
        ).first()
        if existing_frame_check and existing_frame_check.path:
            old_f_dir = Path(existing_frame_check.path).parent
            if old_f_dir.exists() and old_f_dir.resolve() != dest_frame_dir.resolve():
                shutil.rmtree(old_f_dir, ignore_errors=True)

        dest_frame_dir.mkdir(parents=True, exist_ok=True)

        dest_rec_pdb = assert_permanent_path(dest_frame_dir / "receptor.pdb", f"Basin #{rank} Dest Receptor PDB")
        dest_rec_mae = assert_permanent_path(dest_frame_dir / "receptor.mae", f"Basin #{rank} Dest Receptor MAE")
        shutil.copy2(src_rec_pdb, dest_rec_pdb)
        shutil.copy2(src_rec_mae, dest_rec_mae)

        # Idempotent Upsert for ReceptorFrame
        frame_entity = session.query(ReceptorFrame).filter(
            ReceptorFrame.target_id == sys_obj.target_id,
            ReceptorFrame.parent_frame_id == sys_obj.receptor_frame_id,
            ReceptorFrame.feature_group == eff_cluster,
            ReceptorFrame.align_group == eff_align,
            ReceptorFrame.dim_reduction_method == eff_dim,
            ReceptorFrame.cluster_rank == rank,
        ).first()

        if not frame_entity:
            frame_entity = ReceptorFrame(
                target_id=sys_obj.target_id,
                name=new_frame_name,
                source_type="md_centroid",
                parent_frame_id=sys_obj.receptor_frame_id,
                source_trajectory_id=matched_traj.id,
                feature_group=eff_cluster,
                align_group=eff_align,
                dim_reduction_method=eff_dim,
                cluster_rank=rank,
                population=round(float(b_meta["population_pct"]) / 100.0, 4),
                compactness_spread=float(b_meta["compactness_spread"]),
                state="ready",
                path=str(dest_rec_pdb.resolve()),
                path_pdb=str(dest_rec_pdb.resolve()),
                path_mae=str(dest_rec_mae.resolve()),
            )
            session.add(frame_entity)
            session.flush()
        else:
            frame_entity.name = new_frame_name
            frame_entity.source_trajectory_id = matched_traj.id
            frame_entity.align_group = eff_align
            frame_entity.dim_reduction_method = eff_dim
            frame_entity.population = round(float(b_meta["population_pct"]) / 100.0, 4)
            frame_entity.compactness_spread = float(b_meta["compactness_spread"])
            frame_entity.state = "ready"
            frame_entity.path = str(dest_rec_pdb.resolve())
            frame_entity.path_pdb = str(dest_rec_pdb.resolve())
            frame_entity.path_mae = str(dest_rec_mae.resolve())

        centroid_main_pose = None
        centroid_cofactors = []

        for r_code, p_dict in ext_meta.get("poses", {}).items():
            raw_sdf_name = Path(p_dict["sdf_path"]).name
            src_p_sdf = assert_permanent_path(task_dir / raw_sdf_name, f"Basin #{rank} Small Molecule SDF ({r_code})")
            src_p_mae = assert_permanent_path(task_dir / raw_sdf_name.replace(".sdf", ".mae"), f"Basin #{rank} Small Molecule MAE ({r_code})")

            if not src_p_sdf.exists():
                raise FileNotFoundError(f"FATAL: Small molecule pose SDF missing: {src_p_sdf}")

            # Match against database parent by chemical residue name (p_dict["res_name"]), not compound key (r_code)
            base_res_name = p_dict.get("res_name", r_code).strip().upper()
            parent_match = None
            is_main = False
            if sys_obj.pose and base_res_name == sys_obj.pose.res_name.strip().upper():
                parent_match = sys_obj.pose
                is_main = True
            else:
                for cof in sys_obj.cofactors:
                    if base_res_name == cof.res_name.strip().upper():
                        parent_match = cof
                        break

            if not parent_match:
                continue

            sub_inchikey = parent_match.substance.inchikey
            pose_name = f"{new_frame_name}_{r_code}"
            dest_pose_dir = assert_permanent_path(cwd / "substances" / sub_inchikey / "poses" / pose_name, f"Basin #{rank} Pose Directory ({r_code})")
            dest_pose_dir.mkdir(parents=True, exist_ok=True)

            dest_pose_sdf = assert_permanent_path(dest_pose_dir / f"{r_code}_pose.sdf", f"Basin #{rank} Dest Pose SDF ({r_code})")
            dest_pose_mae = assert_permanent_path(dest_pose_dir / f"{r_code}_pose.mae", f"Basin #{rank} Dest Pose MAE ({r_code})")
            shutil.copy2(src_p_sdf, dest_pose_sdf)
            if src_p_mae.exists():
                shutil.copy2(src_p_mae, dest_pose_mae)

            # Idempotent Upsert for Pose
            p_entity = session.query(Pose).filter_by(
                receptor_frame_id=frame_entity.id,
                name=pose_name,
            ).first()

            if not p_entity:
                p_entity = Pose(
                    substance_state_id=parent_match.substance_state_id,
                    receptor_frame_id=frame_entity.id,
                    name=pose_name,
                    res_name=parent_match.res_name,
                    origin="md_cluster",
                    state="ready",
                    path=str(dest_pose_sdf.resolve()),
                    path_sdf=str(dest_pose_sdf.resolve()),
                    path_mae=str(dest_pose_mae.resolve()) if src_p_mae.exists() else None,
                )
                session.add(p_entity)
                session.flush()
            else:
                p_entity.substance_state_id = parent_match.substance_state_id
                p_entity.res_name = parent_match.res_name
                p_entity.path = str(dest_pose_sdf.resolve())
                p_entity.path_sdf = str(dest_pose_sdf.resolve())
                p_entity.path_mae = str(dest_pose_mae.resolve()) if src_p_mae.exists() else None
                p_entity.state = "ready"

            if is_main:
                centroid_main_pose = p_entity
            else:
                centroid_cofactors.append(p_entity)

        # Single canonical 'sys_' prefix applied once to System entity (no 'sys_sys_...')
        sys_centroid_name = f"sys_{new_frame_name}_{sys_obj.pose.substance.name}" if centroid_main_pose else f"sys_{new_frame_name}"
        dest_sys_dir = cwd / "systems" / sys_centroid_name
        dest_sys_dir.mkdir(parents=True, exist_ok=True)

        if src_view_pdb.exists():
            shutil.copy2(src_view_pdb, dest_sys_dir / f"{sys_centroid_name}_view.pdb")

        # Prune any stale systems previously attached to this frame that do not match the new canonical name
        stale_frame_systems = session.query(System).filter(
            System.receptor_frame_id == frame_entity.id,
            System.name != sys_centroid_name,
        ).all()
        for sfs in stale_frame_systems:
            if not sfs.trajectories:
                sfs_dir = cwd / "systems" / sfs.name
                if sfs_dir.exists():
                    shutil.rmtree(sfs_dir, ignore_errors=True)
                session.delete(sfs)

        # Prune any stale poses previously attached to this frame that do not match the new canonical poses
        promoted_pose_names = {f"{new_frame_name}_{r}" for r in ext_meta.get("poses", {}).keys()}
        stale_frame_poses = session.query(Pose).filter(
            Pose.receptor_frame_id == frame_entity.id,
            Pose.name.notin_(promoted_pose_names),
        ).all()
        for sfp in stale_frame_poses:
            if sfp.path:
                sfp_dir = Path(sfp.path).parent
                if sfp_dir.exists() and sfp_dir.name == sfp.name:
                    shutil.rmtree(sfp_dir, ignore_errors=True)
            session.delete(sfp)

        target_sys = session.query(System).filter_by(
            receptor_frame_id=frame_entity.id,
            name=sys_centroid_name,
        ).first()

        if not target_sys:
            target_sys = System(
                target_id=sys_obj.target_id,
                receptor_frame_id=frame_entity.id,
                pose_id=centroid_main_pose.id if centroid_main_pose else None,
                name=sys_centroid_name,
                state="raw",
                origin="md_cluster",
            )
            target_sys.cofactors = centroid_cofactors
            session.add(target_sys)
        else:
            target_sys.pose_id = centroid_main_pose.id if centroid_main_pose else None
            target_sys.cofactors = centroid_cofactors

        promoted += 1

    session.commit()
    session.close()

    console.print(
        f"[bold green]✓ Successfully promoted {promoted} centroid conformation(s) into SQLite 3NF![/bold green]"
    )


def run_cluster(
    systems: list[str] | None = None,
    target: str | None = None,
    traj_names: list[str] | None = None,
    align_group: str | None = None,
    cluster_group: str | None = None,
    dim_reduction: str | None = None,
    n_clusters: int | None = None,
    min_population: float | None = None,
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
    eff_cluster = cluster_group if cluster_group is not None else cluster_cfg["cluster_group"]
    eff_dim = dim_reduction if dim_reduction is not None else cluster_cfg["dim_reduction"]
    eff_k = n_clusters if n_clusters is not None else cluster_cfg["n_clusters"]
    eff_min_pop = min_population if min_population is not None else cluster_cfg["min_cluster_population"]
    eff_auto = auto or cluster_cfg["auto_promote"]

    configs.setdefault("md_simulation", {}).setdefault("clustering", {})["align_group"] = eff_align
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})["cluster_group"] = eff_cluster
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})["dim_reduction"] = eff_dim
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})["n_clusters"] = eff_k
    configs.setdefault("md_simulation", {}).setdefault("clustering", {})["min_cluster_population"] = eff_min_pop

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
        ndx_path = str(Path(sys_obj.path_ndx).resolve()) if sys_obj.path_ndx and Path(sys_obj.path_ndx).exists() else ""
        lig_code = sys_obj.pose.res_name if sys_obj.pose else ""

        if eff_cluster.lower() == "pocket" and not lig_code:
            console.print(
                f"[bold red]FATAL: System '{sys_obj.name}' has no registered ligand pose! "
                "Cannot resolve pocket without reference ligand. Use --cluster-group Backbone or custom residue range.[/bold red]"
            )
            session.close()
            sys.exit(1)

        # Deterministically collect canonical SubstanceState templates for ALL organic components (ligand + cofactors)
        heterogens_spec = []
        template_sdfs = []

        if sys_obj.pose and sys_obj.pose.substance:
            is_ion = (sys_obj.pose.substance.substance_type == "ion")
            p_state = sys_obj.pose.substance_state
            sdf_p = p_state.path_sdf if (p_state and p_state.path_sdf) else sys_obj.pose.path_sdf
            resolved_sdf = str(Path(sdf_p).resolve()) if (sdf_p and Path(sdf_p).exists() and not is_ion) else ""
            if resolved_sdf:
                template_sdfs.append(resolved_sdf)

            heterogens_spec.append({
                "res_name": sys_obj.pose.res_name.strip().upper(),
                "is_ligand": True,
                "is_ion": is_ion,
                "template_sdf": resolved_sdf,
            })

        for cof in sys_obj.cofactors:
            if not cof.substance:
                continue
            is_ion = (cof.substance.substance_type == "ion")
            c_state = cof.substance_state
            sdf_p = c_state.path_sdf if (c_state and c_state.path_sdf) else cof.path_sdf
            resolved_sdf = str(Path(sdf_p).resolve()) if (sdf_p and Path(sdf_p).exists() and not is_ion) else ""
            if resolved_sdf:
                template_sdfs.append(resolved_sdf)

            heterogens_spec.append({
                "res_name": cof.res_name.strip().upper(),
                "is_ligand": False,
                "is_ion": is_ion,
                "template_sdf": resolved_sdf,
            })

        # Purge obsolete centroid deliverables from previous runs with higher K before dispatching
        dest_reports_dir = cwd / "reports" / "clusters" / sys_obj.name / f"{eff_align}_{eff_cluster}"
        if dest_reports_dir.exists():
            for stale_file in dest_reports_dir.glob("centroid_*"):
                stale_file.unlink(missing_ok=True)

        task_id = f"{sys_obj.name}_{eff_align}_{eff_cluster}"

        plan_list.append({
            "task_id": task_id,
            "system": sys_obj.name,
            "clean_xtcs": clean_xtcs,
            "replicate_names": rep_names,
            "ref_pdb": ref_pdb,
            "ndx": ndx_path,
            "align_group": eff_align,
            "cluster_group": eff_cluster,
            "ligand": lig_code,
            "heterogens": heterogens_spec,
            "template_sdfs": list(set(template_sdfs)),
        })

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w", encoding="utf-8") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Dispatching macro-basin clustering (Align: {eff_align}, Features: {eff_cluster}, K={eff_k})...[/bold cyan]"
    )

    receipt_file = cwd / "cluster_receipt.tsv"
    receipt_file.unlink(missing_ok=True)

    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend([
        "--cluster_plan", str(plan_file.resolve()),
    ])

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
                receipt_records.append({
                    "system": parts[0],
                    "align_group": parts[1],
                    "cluster_group": parts[2],
                    "manifest_path": Path(parts[3]).resolve(),
                })

    session.close()

    for rec_meta in receipt_records:
        sys_name = rec_meta["system"]
        m_path = rec_meta["manifest_path"]

        if not m_path.exists() or m_path.stat().st_size == 0:
            raise FileNotFoundError(f"FATAL: Published cluster manifest missing on disk: {m_path}")

        with open(m_path, "r", encoding="utf-8") as jf:
            manifest = json.load(jf)

        basins = manifest.get("basins", [])
        if not basins:
            console.print(f"[yellow]Warning: No basins found in manifest for system '{sys_name}'.[/yellow]")
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
                verdict = "[bold green]Recommended (Major Basin)[/bold green]" if rank == 1 else "[green]Recommended[/green]"
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

        html_report_path = assert_permanent_path(m_path.parent / html_report_name, "Interactive Landscape HTML")
        pml_session_path = assert_permanent_path(m_path.parent / pml_session_name, "Portable PyMOL Session")

        console.print(table)
        console.print(f"📊 [bold cyan]Interactive Landscape (FES):[/bold cyan] {html_report_path}")
        console.print(f"🔍 [bold cyan]Portable PyMOL 3D Session:[/bold cyan] pymol {pml_session_path}\n")

        selected_ranks = []
        if eff_auto or not sys.stdin.isatty():
            selected_ranks = recommended_ranks
            console.print(f"[cyan]Headless Auto-Promote: Selected basins {selected_ranks} for SQLite 3NF registration.[/cyan]")
        else:
            rec_str = ", ".join(map(str, recommended_ranks))
            prompt_msg = f"Action [Enter to accept #{rec_str} | specify ranks (e.g. 1, 3) | (c)ancel]: "
            try:
                user_choice = input(prompt_msg).strip().lower()
            except (EOFError, KeyboardInterrupt):
                console.print("\n[yellow]Clustering review canceled by user.[/yellow]")
                return

            if user_choice in ["c", "cancel"]:
                console.print(
                    f"[yellow]Cluster deliverables preserved in {m_path.parent}. Database not modified.\n"
                    f"To promote later without re-running MD, use: sbdd traj-cluster-promote --manifest {m_path} --ranks 1,2[/yellow]"
                )
                continue
            elif not user_choice:
                selected_ranks = recommended_ranks
            else:
                try:
                    selected_ranks = [int(x.strip().replace("#", "")) for x in user_choice.split(",") if x.strip()]
                except ValueError:
                    console.print("[red]Invalid rank selection. Promotion aborted for this system.[/red]")
                    continue

        if selected_ranks:
            run_cluster_promote(str(m_path), selected_ranks)


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
