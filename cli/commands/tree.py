"""
File: cli/commands/tree.py
Description: CLI handler to visualize the database lineage as a tree structure.
"""

from rich.console import Console
from rich.tree import Tree
from cli.core.db import init_db, Receptor, Ligand, Trajectory

console = Console()


from cli.core.db import init_db, Receptor, Ligand, Pose, System, Trajectory


def cmd_tree(args):
    session = init_db()
    root_tree = Tree("[bold blue]SBDD Project Relational Database[/bold blue]")

    # 1. Systems & Their Trajectories (The Core Execution Units)
    systems = session.query(System).all()
    if systems:
        sys_branch = root_tree.add("[bold green]Systems (Complexes)[/bold green]")
        for s in systems:
            pose_desc = f"[cyan]+ Pose: {s.pose_name}[/cyan]" if s.pose_name else "[dim](Apo)[/dim]"
            s_color = "green" if s.status == "ready" else "yellow"
            s_node = sys_branch.add(f"[{s_color}]{s.name}[/{s_color}] [dim]({s.status})[/dim] : Receptor [cyan]{s.receptor_name}[/cyan] {pose_desc}")

            trajs = session.query(Trajectory).filter_by(system_name=s.name).all()
            for t in trajs:
                t_color = "green" if t.status == "cleaned" else "yellow"
                clean_info = f" -> {t.path_clean_xtc}" if t.path_clean_xtc else f" -> {t.path_xtc}"
                s_node.add(f"[magenta]{t.name}[/magenta] [{t_color}]({t.status})[/{t_color}][dim]{clean_info}[/dim]")

    # 2. Receptors
    receptors = session.query(Receptor).all()
    if receptors:
        rec_branch = root_tree.add("[bold cyan]Receptors[/bold cyan]")
        for r in receptors:
            r_color = "green" if r.status == "ready" else "yellow"
            rec_branch.add(f"[{r_color}]{r.name}[/{r_color}] [dim]({r.status}, {r.prep_mode}) -> {r.path}[/dim]")

    # 3. Ligands & Poses
    ligands = session.query(Ligand).all()
    if ligands:
        lig_branch = root_tree.add("[bold yellow]Ligands & Binding Poses[/bold yellow]")
        for lig in ligands:
            lig_node = lig_branch.add(f"[bold]{lig.name}[/bold] [dim]({lig.status}) -> {lig.path}[/dim]")
            poses = session.query(Pose).filter_by(ligand_name=lig.name).all()
            for p in poses:
                p_color = "green" if p.status == "ready" else "yellow"
                lig_node.add(f"Pose: [{p_color}]{p.name}[/{p_color}] [dim]({p.origin}, in {p.receptor_name}) -> {p.path}[/dim]")

    session.close()
    console.print(root_tree)
