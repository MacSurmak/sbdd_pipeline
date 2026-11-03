"""
File: cli/commands/tree.py
Description: CLI handler to visualize the 3NF database lineage cleanly without string naming hacks.
"""

from pathlib import Path
from rich.console import Console
from rich.tree import Tree
from cli.core.db import (
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
    DockingRun,
    DockedPose,
)

console = Console()


def _format_param_status(st: SubstanceState, sub: Substance) -> str:
    if not (st.path_mol2 and st.path_frcmod):
        return ""
    mol2_p = Path(st.path_mol2)
    frcmod_p = Path(st.path_frcmod)
    frcmod_text = (
        frcmod_p.read_text(encoding="utf-8", errors="ignore")
        if frcmod_p.exists()
        else ""
    )
    is_openff = "OpenFF" in frcmod_text
    opt_sdf = mol2_p.parent / f"{sub.res_name}_opt.sdf"

    if is_openff:
        ff_version = "OpenFF-2.3.0" if "2.3.0" in frcmod_text else "OpenFF"
        if opt_sdf.exists():
            return f" [bold green]+RESP({ff_version})[/bold green]"
        elif "AshGC" in frcmod_text:
            return f" [bold green]+AshGC({ff_version})[/bold green]"
        else:
            return f" [bold green]+BCC({ff_version})[/bold green]"
    else:
        if opt_sdf.exists():
            return " [bold green]+RESP(GAFF2)[/bold green]"
        else:
            return " [bold green]+BCC(GAFF2)[/bold green]"


def _render_detailed_substance_node(parent_node, sub: Substance, session) -> None:
    type_badge = f"[bold]{sub.substance_type.upper()}[/bold]"
    sub_color = "green" if sub.state == "ready" else "yellow"
    states = (
        session.query(SubstanceState)
        .filter_by(substance_id=sub.id)
        .order_by(SubstanceState.state_hash)
        .all()
    )
    smi_suffix = f" [dim]SMILES: {sub.smiles}[/dim]" if not states else ""
    sub_node = parent_node.add(
        f"{type_badge}: [bold cyan]{sub.name}[/bold cyan] [{sub_color}]\\[{sub.state}][/{sub_color}] [dim]({sub.inchikey})[/dim]{smi_suffix}"
    )

    for st in states:
        st_color = "green" if st.state == "ready" else "yellow"
        q_label = (
            f"q={st.formal_charge:+d}" if st.formal_charge != 0 else "neutral"
        )
        param_status = _format_param_status(st, sub)
        state_node = sub_node.add(
            f"State: [magenta]{st.state_hash}[/magenta] [{st_color}]\\[{st.state}][/{st_color}] [cyan]({q_label})[/cyan]{param_status} [dim]SMILES: {st.isomeric_smiles}[/dim]"
        )

        poses = (
            session.query(Pose)
            .filter_by(substance_state_id=st.id)
            .order_by(Pose.name)
            .all()
        )
        for p in poses:
            p_color = "green" if p.state == "ready" else "yellow"
            state_node.add(
                f"Pose: [{p_color}]{p.name}[/{p_color}] [dim]({p.origin}, state={p.state})[/dim] in [cyan]{p.receptor_frame.name}[/cyan] -> {Path(p.path).name}"
            )


def print_tree(
    batch: str | None = None,
    substance: str | None = None,
) -> None:
    from sqlalchemy import func

    session = init_db()

    # Mode A: Targeted single-substance inspection
    if substance:
        query_val = substance.strip()
        sub_obj = (
            session.query(Substance)
            .filter((Substance.name == query_val) | (Substance.inchikey == query_val))
            .first()
        )
        if not sub_obj:
            console.print(
                f"[bold red]FATAL: Substance '{query_val}' not found in database.[/bold red]"
            )
            session.close()
            return

        sub_root = Tree(
            f"[bold blue]Substance Lineage Inspection: {sub_obj.name}[/bold blue] [dim]({sub_obj.inchikey})[/dim]"
        )
        sub_root.add(f"Batch Provenance: [bold cyan]{sub_obj.batch_id or 'none'}[/bold cyan]")
        _render_detailed_substance_node(sub_root, sub_obj, session)

        # Show evaluation across DockingRuns
        d_poses = (
            session.query(DockedPose)
            .filter_by(substance_id=sub_obj.id)
            .order_by(DockedPose.docking_run_id, DockedPose.mode_index)
            .all()
        )
        if d_poses:
            dock_branch = sub_root.add("[bold magenta]Analytical Docking Evaluations[/bold magenta]")
            by_run = {}
            for dp in d_poses:
                by_run.setdefault((dp.docking_run.name, dp.receptor_frame.name), []).append(dp)
            for (r_name, f_name), modes in by_run.items():
                g_modes = [m for m in modes if m.engine == "glide"]
                n_modes = [m for m in modes if m.engine == "gnina"]
                g_best = f"GScore: {g_modes[0].glide_score:.2f}" if g_modes and g_modes[0].glide_score is not None else ""
                n_best = f"CNNaff: {n_modes[0].cnn_affinity:.2f}" if n_modes and n_modes[0].cnn_affinity is not None else ""
                stats_str = ", ".join(s for s in [g_best, n_best] if s)
                stats_label = f" [dim]({stats_str})[/dim]" if stats_str else ""
                dock_branch.add(
                    f"Run [cyan]{r_name}[/cyan] on Frame [cyan]{f_name}[/cyan]: {len(modes)} mode(s) (Glide: {len(g_modes)}, GNINA: {len(n_modes)}){stats_label}"
                )

        session.close()
        console.print(sub_root)
        return

    root_tree = Tree("[bold blue]SBDD Project Relational Database (3NF)[/bold blue]")

    # 1. Biological Targets, Frames & Systems
    targets = session.query(Target).order_by(Target.name).all()
    if targets:
        t_branch = root_tree.add(
            "[bold green]Biological Targets & Systems[/bold green]"
        )
        for target in targets:
            target_node = t_branch.add(
                f"[bold magenta]Target:[/bold magenta] {target.name}"
            )

            # Frames
            frames = (
                session.query(ReceptorFrame)
                .filter_by(target_id=target.id)
                .order_by(ReceptorFrame.name)
                .all()
            )
            if frames:
                f_branch = target_node.add(
                    "[bold cyan]Receptor Frames (3D)[/bold cyan]"
                )
                for f in frames:
                    f_color = "green" if f.state == "ready" else "yellow"
                    lineage = (
                        f" [dim](derived from #{f.parent_frame_id})[/dim]"
                        if f.parent_frame_id
                        else ""
                    )
                    f_branch.add(
                        f"[{f_color}]{f.name}[/{f_color}] [dim]({f.source_type}, state={f.state})[/dim]{lineage}"
                    )

            # Systems
            systems = (
                session.query(System)
                .filter_by(target_id=target.id)
                .order_by(System.name)
                .all()
            )
            if systems:
                s_branch = target_node.add("[bold blue]Physical Systems[/bold blue]")
                for s in systems:
                    s_color = "green" if s.state == "ready" else "yellow"
                    pose_label = (
                        f" + Ligand [cyan]{s.pose.name}[/cyan]"
                        if s.pose
                        else " [dim](Apo)[/dim]"
                    )
                    cofs = [c.substance.name for c in s.cofactors if c.substance]
                    cofs_label = (
                        f" [yellow]+[{','.join(cofs)}][/yellow]" if cofs else ""
                    )
                    sys_node = s_branch.add(
                        f"[{s_color}]{s.name}[/{s_color}] [dim](state={s.state})[/dim] on [cyan]{s.receptor_frame.name}[/cyan]{pose_label}{cofs_label}"
                    )

                    # Trajectories under system
                    trajs = (
                        session.query(Trajectory)
                        .filter_by(system_id=s.id)
                        .order_by(Trajectory.replicate)
                        .all()
                    )
                    if trajs:
                        tr_branch = sys_node.add("[dim]Trajectories[/dim]")
                        for t in trajs:
                            if t.status == "cleaned":
                                t_color = "green"
                            elif t.status in ["drifted", "failed", "disabled"]:
                                t_color = "bold red"
                            else:
                                t_color = "yellow"
                            tr_branch.add(
                                f"[magenta]{t.name}[/magenta] [{t_color}]({t.status})[/{t_color}] - {t.length_ns} ns"
                            )

                    # MM-PBSA Runs under system
                    mmpbsa_runs = (
                        session.query(MmpbsaRun)
                        .filter_by(system_id=s.id)
                        .order_by(MmpbsaRun.created_at)
                        .all()
                    )
                    if mmpbsa_runs:
                        mmpbsa_branch = sys_node.add(
                            "[bold magenta]MM-PBSA Thermodynamic Runs[/bold magenta]"
                        )
                        for m_run in mmpbsa_runs:
                            boltz_str = (
                                f", Boltzmann ΔG={m_run.delta_g_boltzmann_total:.2f} kcal/mol"
                                if m_run.delta_g_boltzmann_total is not None
                                else ""
                            )
                            mmpbsa_branch.add(
                                f"[magenta]Run #{m_run.id}[/magenta]: [green]ΔG={m_run.delta_g_overall_mean:.2f} ± {m_run.delta_g_overall_sem:.2f} kcal/mol[/green]"
                                f"{boltz_str} [dim](interval={m_run.sample_interval_ps} ps, T={m_run.temperature} K)[/dim]"
                            )

                    # Clustering Runs under system
                    c_runs = (
                        session.query(ClusteringRun)
                        .filter_by(system_id=s.id)
                        .order_by(ClusteringRun.created_at)
                        .all()
                    )
                    if c_runs:
                        c_branch = sys_node.add(
                            "[bold cyan]Conformational Clusterings[/bold cyan]"
                        )
                        for cr in c_runs:
                            cr_node = c_branch.add(
                                f"[cyan]{cr.name}[/cyan] [dim]({cr.algorithm}, {cr.feature_group} on {cr.align_group}, {len(cr.clusters)} basins)[/dim]"
                            )
                            for cl in cr.clusters:
                                dg_str = (
                                    f" -> [bold green]ΔG={cl.delta_g_mean:.2f} ± {cl.delta_g_sem:.2f} kcal/mol[/bold green]"
                                    if cl.delta_g_mean is not None
                                    else ""
                                )
                                cr_node.add(
                                    f"Basin #{cl.cluster_rank} ({cl.population_fraction * 100:.1f}%, {cl.frame_count} frames, medoid: {cl.medoid_replicate}@{cl.medoid_time_ns:.1f}ns){dg_str}"
                                )

    # 1.5 Analytical Docking Runs
    dock_runs = session.query(DockingRun).order_by(DockingRun.created_at).all()
    if dock_runs:
        d_branch = root_tree.add(
            "[bold magenta]Analytical Docking Runs (3NF)[/bold magenta]"
        )
        for dr in dock_runs:
            total_poses = len(dr.docked_poses)
            frames_set = sorted(
                list({dp.receptor_frame.name for dp in dr.docked_poses})
            )
            subs_set = sorted(list({dp.substance.name for dp in dr.docked_poses}))
            dr_node = d_branch.add(
                f"[bold cyan]{dr.name}[/bold cyan] [dim](mode={dr.mode}, {len(frames_set)} frame(s) × {len(subs_set)} substance(s), {total_poses} total modes)[/dim]"
            )
            for f_name in frames_set:
                f_poses = [
                    dp for dp in dr.docked_poses if dp.receptor_frame.name == f_name
                ]
                g_count = sum(1 for dp in f_poses if dp.engine == "glide")
                n_count = sum(1 for dp in f_poses if dp.engine == "gnina")
                dr_node.add(
                    f"Frame [cyan]{f_name}[/cyan]: {len( {dp.substance_id for dp in f_poses} )} substance(s) [dim](Glide: {g_count} modes, GNINA: {n_count} modes)[/dim]"
                )

    # 2. Chemical Substances: Structural Focus Entities vs Screening Library Pools
    sub_branch = root_tree.add(
        "[bold yellow]Chemical Substances & Microstates (3NF)[/bold yellow]"
    )

    # 2A. Structural and Active Substances (always rendered in full detail)
    # Includes crystal/medoid extracts and any substance that has physical Poses
    structural_substances = (
        session.query(Substance)
        .filter(
            (Substance.batch_id == "crystal_extract")
            | (
                Substance.id.in_(
                    session.query(SubstanceState.substance_id)
                    .join(Pose, Pose.substance_state_id == SubstanceState.id)
                )
            )
        )
        .order_by(Substance.name)
        .all()
    )
    if structural_substances:
        struct_branch = sub_branch.add(
            f"[bold green]Structural & Bound Entities[/bold green] [dim]({len(structural_substances)} entity(ies))[/dim]"
        )
        for sub in structural_substances:
            _render_detailed_substance_node(struct_branch, sub, session)

    # 2B. Screening Library Pools (batch_id != 'crystal_extract')
    # Filter by --batch if requested
    batch_filter = batch.strip() if batch else None

    # Fast SQL aggregation per batch: total substances, ready count, raw count
    stats_query = (
        session.query(
            Substance.batch_id,
            func.count(Substance.id).label("total"),
            func.sum(func.case((Substance.state == "ready", 1), else_=0)).label("ready_cnt"),
            func.sum(func.case((Substance.state == "raw", 1), else_=0)).label("raw_cnt"),
        )
        .filter(Substance.batch_id != "crystal_extract")
        .group_by(Substance.batch_id)
    )
    if batch_filter:
        stats_query = stats_query.filter(Substance.batch_id == batch_filter)

    batch_stats = stats_query.order_by(Substance.batch_id).all()

    if batch_stats:
        pool_branch = sub_branch.add(
            f"[bold cyan]Screening Library Pools[/bold cyan] [dim]({len(batch_stats)} batch(es))[/dim]"
        )

        for b_name, total_subs, ready_subs, raw_subs in batch_stats:
            b_tag = b_name or "unspecified_batch"
            all_ready = (ready_subs == total_subs)
            b_style = "green" if all_ready else "yellow"
            b_status_str = f"[{b_style}]\\[{'ready' if all_ready else 'raw'}][/{b_style}]"

            # Query docking evaluations for this batch
            docking_summary = (
                session.query(
                    DockingRun.name,
                    func.count(func.distinct(DockedPose.substance_id)),
                    func.count(DockedPose.id),
                )
                .join(DockedPose, DockedPose.docking_run_id == DockingRun.id)
                .join(Substance, DockedPose.substance_id == Substance.id)
                .filter(Substance.batch_id == b_name)
                .group_by(DockingRun.name)
                .all()
            )

            # If user explicitly requested --batch <name> or batch is tiny (<= 5), expand substances
            should_expand = bool(batch_filter) or (total_subs <= 5)

            batch_node = pool_branch.add(
                f"[bold cyan]Batch: {b_tag}[/bold cyan] {b_status_str} [dim]({total_subs} substance(s): {ready_subs} ready, {raw_subs} raw)[/dim]"
            )

            if docking_summary:
                d_info = ", ".join(
                    f"'{run_n}' ({s_cnt}/{total_subs} subs, {m_cnt} modes)"
                    for run_n, s_cnt, m_cnt in docking_summary
                )
                batch_node.add(f"[magenta]Docking Evaluation:[/magenta] {d_info}")
            else:
                batch_node.add("[dim]Docking Evaluation: not evaluated[/dim]")

            if should_expand:
                subs_in_batch = (
                    session.query(Substance)
                    .filter_by(batch_id=b_name)
                    .order_by(Substance.name)
                    .all()
                )
                for sub in subs_in_batch:
                    _render_detailed_substance_node(batch_node, sub, session)
            else:
                batch_node.add(
                    f"[dim]Run 'sbdd tree --batch {b_tag}' to inspect all {total_subs} substances or 'sbdd tree --substance <name>' for a specific molecule.[/dim]"
                )

    session.close()
    console.print(root_tree)

    session.close()
    console.print(root_tree)
