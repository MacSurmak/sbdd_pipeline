"""
File: cli/commands/tree.py
Description: CLI handler to visualize the 3NF database lineage cleanly without string naming hacks.
"""

from pathlib import Path
from rich.console import Console
from rich.tree import Tree
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System, Trajectory

console = Console()


def cmd_tree(args):
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
                    cofs = (
                        f" [yellow]+{','.join(f.cofactors_list)}[/yellow]"
                        if f.cofactors_list
                        else ""
                    )
                    f_branch.add(
                        f"[{f_color}]{f.name}[/{f_color}] [dim]({f.source_type}, state={f.state})[/dim]{cofs}{lineage}"
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
                        f" + Pose [cyan]{s.pose.name}[/cyan]"
                        if s.pose
                        else " [dim](Apo)[/dim]"
                    )
                    sys_node = s_branch.add(
                        f"[{s_color}]{s.name}[/{s_color}] [dim](state={s.state})[/dim] on [cyan]{s.receptor_frame.name}[/cyan]{pose_label}"
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

    # 2. Chemical Substances & Their 3D Poses
    ligands = session.query(Ligand).order_by(Ligand.name).all()
    if ligands:
        lig_branch = root_tree.add(
            "[bold yellow]Chemical Ligands & 3D Poses[/bold yellow]"
        )
        for lig in ligands:
            lig_node = lig_branch.add(
                f"[bold yellow]{lig.name}[/bold yellow] [dim]({lig.inchikey}) SMILES: {lig.smiles}[/dim]"
            )
            poses = (
                session.query(Pose)
                .filter_by(ligand_id=lig.id)
                .order_by(Pose.name)
                .all()
            )
            for p in poses:
                p_color = "green" if p.state == "ready" else "yellow"
                ch = p.scores_dict.get("formal_charge", None)
                ch_tag = (
                    f" [cyan](q={ch:+d})[/cyan]"
                    if ch is not None and ch != 0
                    else (" [dim](neutral)[/dim]" if ch == 0 else "")
                )
                lig_node.add(
                    f"Pose: [{p_color}]{p.name}[/{p_color}] [dim]({p.origin}, state={p.state})[/dim]{ch_tag} in [cyan]{p.receptor_frame.name}[/cyan] -> {Path(p.path).name}"
                )

    session.close()
    console.print(root_tree)
