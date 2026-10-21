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
)

console = Console()


def print_tree() -> None:
    session = init_db()
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

    # 2. Chemical Substances, Microstates & Their 3D Poses
    substances = (
        session.query(Substance)
        .order_by(Substance.substance_type, Substance.name)
        .all()
    )
    if substances:
        sub_branch = root_tree.add(
            "[bold yellow]Chemical Substances & Microstates (3NF)[/bold yellow]"
        )
        for sub in substances:
            type_badge = f"[bold]{sub.substance_type.upper()}[/bold]"
            sub_node = sub_branch.add(
                f"{type_badge}: [bold cyan]{sub.name}[/bold cyan] [dim]({sub.inchikey})[/dim]"
            )

            states = (
                session.query(SubstanceState)
                .filter_by(substance_id=sub.id)
                .order_by(SubstanceState.state_hash)
                .all()
            )
            for st in states:
                st_color = "green" if st.state == "ready" else "yellow"
                q_label = (
                    f"q={st.formal_charge:+d}" if st.formal_charge != 0 else "neutral"
                )
                param_status = ""
                if st.path_mol2 and st.path_frcmod:
                    # Detect if ORCA RESP was used by inspecting relaxed geometry artifact
                    mol2_p = Path(st.path_mol2)
                    opt_sdf = mol2_p.parent / f"{sub.res_name}_opt.sdf"
                    if opt_sdf.exists():
                        param_status = " [bold green]+RESP(GAFF2)[/bold green]"
                    else:
                        param_status = " [bold green]+BCC(GAFF2)[/bold green]"

                state_node = sub_node.add(
                    f"State: [magenta]{st.state_hash}[/magenta] [{st_color}][{st.state}][/{st_color}] [cyan]({q_label})[/cyan]{param_status} [dim]SMILES: {st.isomeric_smiles}[/dim]"
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

    session.close()
    console.print(root_tree)
