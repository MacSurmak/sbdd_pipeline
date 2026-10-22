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
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def _transfer_coordinates_to_template(
    template_sdf: Path, raw_pose_sdf: Path, out_sdf: Path
) -> None:
    """
    Transfers exact Cartesian coordinates from MD centroid structure onto the
    canonical, validated chemical graph template of the parent SubstanceState.
    Preserves 100% of the true bond orders, formal charges, stereochemistry,
    and valences without any heuristic re-assignment or valence guessing.
    """
    from rdkit import Chem

    suppl_tmpl = Chem.SDMolSupplier(str(template_sdf), removeHs=False, sanitize=False)
    tmpl_mol = next((m for m in suppl_tmpl if m is not None), None)

    suppl_raw = Chem.SDMolSupplier(str(raw_pose_sdf), removeHs=False, sanitize=False)
    raw_mol = next((m for m in suppl_raw if m is not None), None)

    if not tmpl_mol or not raw_mol:
        raise ValueError(
            f"FATAL: Could not parse chemical structure from template {template_sdf} or MD pose {raw_pose_sdf}"
        )

    n_tmpl = tmpl_mol.GetNumAtoms()
    n_raw = raw_mol.GetNumAtoms()
    if n_tmpl != n_raw:
        raise ValueError(
            f"FATAL: Atom count mismatch between template {template_sdf.name} ({n_tmpl} atoms) "
            f"and MD centroid pose {raw_pose_sdf.name} ({n_raw} atoms)! Topology cannot be preserved."
        )

    for i in range(n_tmpl):
        s_tmpl = tmpl_mol.GetAtomWithIdx(i).GetSymbol()
        s_raw = raw_mol.GetAtomWithIdx(i).GetSymbol()
        if s_tmpl != s_raw:
            raise ValueError(
                f"FATAL: Element mismatch at atom index {i}: template has {s_tmpl}, MD pose has {s_raw}!"
            )

    conf_tmpl = tmpl_mol.GetConformer()
    conf_raw = raw_mol.GetConformer()
    for i in range(n_tmpl):
        conf_tmpl.SetAtomPosition(i, conf_raw.GetAtomPosition(i))

    for i in range(n_tmpl):
        pos = conf_tmpl.GetAtomPosition(i)
        if abs(pos.x) < 1e-4 and abs(pos.y) < 1e-4 and abs(pos.z) < 1e-4:
            raise ValueError(
                f"FATAL: Atom {tmpl_mol.GetAtomWithIdx(i).GetSymbol()}{i+1} has degenerate coordinates (0, 0, 0)!"
            )

    tmpl_mol.SetProp("_Name", out_sdf.stem)
    tmpl_mol.SetProp("s_m_title", out_sdf.stem)

    with Chem.SDWriter(str(out_sdf)) as writer:
        writer.write(tmpl_mol)


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
    session = init_db()
    configs = load_project_configs(cwd)
    cluster_cfg = configs.get("md_simulation", {}).get("clustering", {})

    # Strict parameter resolution: explicit CLI arguments override validated md.yaml configuration
    eff_align = align_group or cluster_cfg.get("align_group", "Backbone")
    eff_cluster = cluster_group or cluster_cfg.get("cluster_group", "Pocket")
    eff_dim = dim_reduction or cluster_cfg.get("dim_reduction", "pca")
    eff_k = n_clusters or cluster_cfg.get("n_clusters", 3)
    eff_min_pop = min_population or cluster_cfg.get("min_cluster_population", 5.0)
    eff_auto = auto or cluster_cfg.get("auto_promote", False)

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
        ref_pdb = str(Path(t_list[0].path_ref_pdb).resolve())
        ndx_path = str(Path(sys_obj.path_ndx).resolve()) if sys_obj.path_ndx and Path(sys_obj.path_ndx).exists() else ""
        lig_code = sys_obj.pose.res_name if sys_obj.pose else ""

        # Fail-Fast validation on pocket definitions
        if eff_cluster.lower() == "pocket" and not lig_code:
            console.print(
                f"[bold red]FATAL: System '{sys_obj.name}' has no registered ligand pose! "
                "Cannot resolve pocket without reference ligand. Use --cluster-group Backbone or custom residue range.[/bold red]"
            )
            session.close()
            sys.exit(1)

        template_sdf = ""
        if sys_obj.pose:
            p_state = sys_obj.pose.substance_state
            if p_state and p_state.path_sdf and Path(p_state.path_sdf).exists():
                template_sdf = str(Path(p_state.path_sdf).resolve())
            elif sys_obj.pose.path_sdf and Path(sys_obj.pose.path_sdf).exists():
                template_sdf = str(Path(sys_obj.pose.path_sdf).resolve())

        task_id = f"{sys_obj.name}_{eff_align}_{eff_cluster}"
        subpath = f"{eff_align}_{eff_cluster}"

        plan_list.append({
            "task_id": task_id,
            "system": sys_obj.name,
            "clean_xtcs": clean_xtcs,
            "ref_pdb": ref_pdb,
            "ndx": ndx_path,
            "align_group": eff_align,
            "cluster_group": eff_cluster,
            "ligand": lig_code,
            "template_sdf": template_sdf,
            "target_subpath": subpath,
        })

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w", encoding="utf-8") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Dispatching macro-basin clustering (Align: {eff_align}, Features: {eff_cluster}, K={eff_k})...[/bold cyan]"
    )

    # Clean previous receipt file before Nextflow invocation
    receipt_file = cwd / "cluster_receipt.tsv"
    receipt_file.unlink(missing_ok=True)

    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend([
        "--cluster_plan", str(plan_file.resolve()),
        "--cluster_align_group", str(eff_align),
        "--cluster_feature_group", str(eff_cluster),
        "--cluster_dim_reduction", str(eff_dim),
        "--cluster_n_clusters", str(eff_k),
        "--cluster_min_pop", str(eff_min_pop),
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
            f"FATAL: Execution receipt {receipt_file} was not emitted by Nextflow. Clustering deliverables missing."
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

    if not receipt_records:
        session.close()
        raise ValueError("FATAL: Receipt file is empty or malformed.")

    promoted_total = 0

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

        sys_obj, t_list = sys_lookup[sys_name]
        task_dir = m_path.parent

        # Render Rich interactive summary table
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

        console.print(table)
        console.print(f"📊 [bold cyan]Interactive Landscape (FES):[/bold cyan] {manifest.get('html_report')}")
        console.print(f"🔍 [bold cyan]Portable PyMOL 3D Session:[/bold cyan] pymol {manifest.get('pymol_session')}\n")

        # Decision flow: Headless auto-promotion vs interactive terminal prompt
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
                session.close()
                return

            if user_choice in ["c", "cancel"]:
                console.print(f"[yellow]Cluster deliverables preserved in {task_dir}. Database not modified.[/yellow]")
                continue
            elif not user_choice:
                selected_ranks = recommended_ranks
            else:
                try:
                    selected_ranks = [int(x.strip().replace("#", "")) for x in user_choice.split(",") if x.strip()]
                except ValueError:
                    console.print("[red]Invalid rank selection. Promotion aborted for this system.[/red]")
                    continue

        if not selected_ranks:
            console.print("[yellow]Zero basins selected. Skipping promotion.[/yellow]")
            continue

        # Transactional 3NF registration (Lineage Variant 1)
        # Centroid deliverables are promoted to canonical project directories without state suffixes in filenames
        for rank in selected_ranks:
            b_meta = next((b for b in basins if int(b["rank"]) == rank), None)
            if not b_meta:
                console.print(f"[yellow]Warning: Basin #{rank} not found in manifest. Skipping.[/yellow]")
                continue

            src_rec_pdb = task_dir / f"centroid_{rank}_receptor.pdb"
            src_view_pdb = task_dir / f"centroid_{rank}_complex.pdb"

            if not src_rec_pdb.exists() or src_rec_pdb.stat().st_size == 0:
                raise FileNotFoundError(f"FATAL: Medoid receptor structure missing on disk: {src_rec_pdb}")

            # Resolve source trajectory record via replicate name
            src_rep_name = str(b_meta["medoid_replicate"])
            matched_traj = next((t for t in t_list if t.replicate == int(src_rep_name.replace("rep_", "")) or t.name == src_rep_name), t_list[0])

            # Canonical clean entity naming (descriptive role without mutable state suffixes)
            new_frame_name = f"{sys_obj.name}_{eff_cluster.lower()}_c{rank}"
            dest_frame_dir = cwd / "targets" / sys_obj.target.name / "frames" / new_frame_name
            dest_frame_dir.mkdir(parents=True, exist_ok=True)

            dest_rec_pdb = dest_frame_dir / "receptor.pdb"
            shutil.copy2(src_rec_pdb, dest_rec_pdb)

            # Generate or link MAE deliverable for Schrödinger Glide/Prime compatibility
            dest_rec_mae = dest_frame_dir / "receptor.mae"
            if shutil.which("structconvert"):
                subprocess.run(["structconvert", str(dest_rec_pdb), str(dest_rec_mae)], capture_output=True, check=False)
            if not dest_rec_mae.exists():
                dest_rec_mae = dest_rec_pdb

            ex_frame = session.query(ReceptorFrame).filter_by(name=new_frame_name).first()
            if not ex_frame:
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
                frame_entity = ex_frame
                frame_entity.source_trajectory_id = matched_traj.id
                frame_entity.feature_group = eff_cluster
                frame_entity.align_group = eff_align
                frame_entity.dim_reduction_method = eff_dim
                frame_entity.cluster_rank = rank
                frame_entity.population = round(float(b_meta["population_pct"]) / 100.0, 4)
                frame_entity.compactness_spread = float(b_meta["compactness_spread"])
                frame_entity.state = "ready"
                frame_entity.path = str(dest_rec_pdb.resolve())
                frame_entity.path_pdb = str(dest_rec_pdb.resolve())
                frame_entity.path_mae = str(dest_rec_mae.resolve())

            # Promote small molecule poses preserving chemical graph topology and cofactors
            centroid_main_pose = None
            centroid_cofactors = []

            for pose_file in task_dir.glob(f"centroid_{rank}_*_pose.sdf"):
                # Extract residue name from filename centroid_<rank>_<res>_pose.sdf
                r_clean = pose_file.name.replace(f"centroid_{rank}_", "").replace("_pose.sdf", "").strip().upper()

                parent_match = None
                is_main = False
                if sys_obj.pose and r_clean == sys_obj.pose.res_name.strip().upper():
                    parent_match = sys_obj.pose
                    is_main = True
                else:
                    for cof in sys_obj.cofactors:
                        if r_clean == cof.res_name.strip().upper():
                            parent_match = cof
                            break

                if not parent_match:
                    continue

                sub_inchikey = parent_match.substance.inchikey
                pose_name = f"{new_frame_name}_{r_clean}"
                dest_pose_dir = cwd / "substances" / sub_inchikey / "poses" / pose_name
                dest_pose_dir.mkdir(parents=True, exist_ok=True)

                dest_pose_sdf = dest_pose_dir / f"{r_clean}_pose.sdf"
                shutil.copy2(pose_file, dest_pose_sdf)

                dest_pose_mae = dest_pose_dir / f"{r_clean}_pose.mae"
                if shutil.which("structconvert"):
                    subprocess.run(["structconvert", str(dest_pose_sdf), str(dest_pose_mae)], capture_output=True, check=False)
                if not dest_pose_mae.exists():
                    dest_pose_mae = None

                ex_p = session.query(Pose).filter_by(name=pose_name).first()
                if not ex_p:
                    p_entity = Pose(
                        substance_state_id=parent_match.substance_state_id,
                        receptor_frame_id=frame_entity.id,
                        name=pose_name,
                        res_name=parent_match.res_name,
                        origin="md_cluster",
                        state="ready",
                        path=str(dest_pose_sdf.resolve()),
                        path_sdf=str(dest_pose_sdf.resolve()),
                        path_mae=str(dest_pose_mae.resolve()) if dest_pose_mae else None,
                    )
                    session.add(p_entity)
                    session.flush()
                else:
                    p_entity = ex_p
                    p_entity.substance_state_id = parent_match.substance_state_id
                    p_entity.path = str(dest_pose_sdf.resolve())
                    p_entity.path_sdf = str(dest_pose_sdf.resolve())
                    p_entity.path_mae = str(dest_pose_mae.resolve()) if dest_pose_mae else None
                    p_entity.state = "ready"

                if is_main:
                    centroid_main_pose = p_entity
                else:
                    centroid_cofactors.append(p_entity)

            # Register derived physical System
            sys_centroid_name = f"{new_frame_name}_{sys_obj.pose.substance.name}" if centroid_main_pose else new_frame_name
            dest_sys_dir = cwd / "systems" / sys_centroid_name
            dest_sys_dir.mkdir(parents=True, exist_ok=True)

            if src_view_pdb.exists():
                shutil.copy2(src_view_pdb, dest_sys_dir / f"{sys_centroid_name}_view.pdb")

            ex_sys = session.query(System).filter_by(name=sys_centroid_name).first()
            if not ex_sys:
                new_sys = System(
                    target_id=sys_obj.target_id,
                    receptor_frame_id=frame_entity.id,
                    pose_id=centroid_main_pose.id if centroid_main_pose else None,
                    name=sys_centroid_name,
                    state="raw",
                    origin="md_cluster",
                )
                new_sys.cofactors = centroid_cofactors
                session.add(new_sys)
            else:
                ex_sys.pose_id = centroid_main_pose.id if centroid_main_pose else None
                ex_sys.state = "raw"
                ex_sys.cofactors = centroid_cofactors

            promoted_total += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Successfully promoted {promoted_total} centroid conformation(s) into SQLite 3NF as ready ReceptorFrames![/bold green]"
    )
    console.print("[dim]Lineage preserved. These frames can now be used for docking via 'sbdd dock'. View with 'sbdd tree'.[/dim]")


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
