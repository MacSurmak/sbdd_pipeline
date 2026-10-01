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

    # 1. Systems & Their Conformational Ensembles
    systems = session.query(System).all()
    if systems:
        sys_branch = root_tree.add(
            "[bold green]Systems (Complexes & Ensembles)[/bold green]"
        )
        root_systems = [s for s in systems if s.parent_system_name is None]

        def render_system(node, sys_obj):
            s_color = "green" if sys_obj.status == "ready" else "yellow"
            pose_info = (
                f"[cyan]+ Pose: {sys_obj.pose_name}[/cyan]"
                if sys_obj.pose_name
                else "[dim](Apo)[/dim]"
            )
            sys_node = node.add(
                f"[{s_color}]{sys_obj.name}[/{s_color}] [dim]({sys_obj.status}, {sys_obj.origin})[/dim] : Receptor [cyan]{sys_obj.receptor_name}[/cyan] {pose_info}"
            )

            # Direct trajectories of this system
            trajs = session.query(Trajectory).filter_by(system_name=sys_obj.name).all()
            if trajs:
                traj_branch = sys_node.add("[dim]Trajectories[/dim]")
                for t in trajs:
                    t_color = "green" if t.status == "cleaned" else "yellow"
                    clean_info = (
                        f" -> {t.path_clean_xtc}"
                        if t.path_clean_xtc
                        else f" -> {t.path_xtc}"
                    )
                    traj_branch.add(
                        f"[magenta]{t.name}[/magenta] [{t_color}]({t.status})[/{t_color}][dim]{clean_info}[/dim]"
                    )

            # Alternative conformations (MD cluster centroids)
            child_confs = [s for s in systems if s.parent_system_name == sys_obj.name]
            if child_confs:
                conf_branch = sys_node.add(
                    "[bold cyan]Alternative Conformations (MD Centroids)[/bold cyan]"
                )
                for child_sys in child_confs:
                    render_system(conf_branch, child_sys)

        for r_sys in root_systems:
            render_system(sys_branch, r_sys)

    # 2. Receptors
    receptors = session.query(Receptor).all()
    if receptors:
        rec_branch = root_tree.add(
            "[bold cyan]Receptors (Conformational Models)[/bold cyan]"
        )
        for r in receptors:
            r_color = "green" if r.status == "ready" else "yellow"
            parent_info = (
                f" [dim](parent: {r.parent_name})[/dim]" if r.parent_name else ""
            )
            rec_branch.add(
                f"[{r_color}]{r.name}[/{r_color}] [dim]({r.status}, {r.prep_mode}) -> {r.path}[/dim]{parent_info}"
            )

    # 3. Ligands & Binding Poses
    ligands = session.query(Ligand).all()
    if ligands:
        lig_branch = root_tree.add("[bold yellow]Ligands & Binding Poses[/bold yellow]")
        for lig in ligands:
            lig_node = lig_branch.add(
                f"[bold]{lig.name}[/bold] [dim]({lig.status}) -> {lig.path}[/dim]"
            )
            poses = session.query(Pose).filter_by(ligand_name=lig.name).all()
            for p in poses:
                p_color = "green" if p.status == "ready" else "yellow"
                lig_node.add(
                    f"Pose: [{p_color}]{p.name}[/{p_color}] [dim]({p.origin}, in {p.receptor_name}) -> {p.path}[/dim]"
                )

    session.close()
    console.print(root_tree)
