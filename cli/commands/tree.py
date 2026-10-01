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

    def get_root_target(rec_name):
        rec = session.query(Receptor).filter_by(name=rec_name).first()
        while rec and rec.parent_name:
            rec = session.query(Receptor).filter_by(name=rec.parent_name).first()
        return rec.name if rec else rec_name

    # 1. Biological Grouping: Target -> Ligand -> Systems (Conformations)
    systems = session.query(System).all()
    if systems:
        tree_data = {}
        for s in systems:
            target = get_root_target(s.receptor_name)
            lig = "Apo"
            if s.pose_name:
                pose = session.query(Pose).filter_by(name=s.pose_name).first()
                if pose:
                    lig = pose.ligand_name
            tree_data.setdefault(target, {}).setdefault(lig, []).append(s)

        sys_branch = root_tree.add("[bold green]Biological Targets & Conformations[/bold green]")
        
        for target, ligs in tree_data.items():
            target_node = sys_branch.add(f"[bold magenta]Target:[/bold magenta] {target}")
            
            for lig, sys_list in ligs.items():
                lig_color = "cyan" if lig != "Apo" else "dim"
                lig_node = target_node.add(f"[{lig_color}]Ligand: {lig}[/{lig_color}]")
                
                for sys_obj in sys_list:
                    s_color = "green" if sys_obj.status == "ready" else "yellow"
                    pose_info = f"[dim](Pose: {sys_obj.pose_name})[/dim]" if sys_obj.pose_name else ""
                    sys_node = lig_node.add(
                        f"[{s_color}]{sys_obj.name}[/{s_color}] [dim]({sys_obj.origin})[/dim] {pose_info}"
                    )
                    
                    trajs = session.query(Trajectory).filter_by(system_name=sys_obj.name).all()
                    if trajs:
                        traj_branch = sys_node.add("[dim]Trajectories[/dim]")
                        for t in trajs:
                            t_color = "green" if t.status == "cleaned" else "yellow"
                            len_info = f" ({t.length_ns} ns)" if t.length_ns else ""
                            traj_branch.add(f"[magenta]{t.name}[/magenta] [{t_color}]({t.status})[/{t_color}]{len_info}")

    # 2. Receptors (Raw structures & intermediates)
    receptors = session.query(Receptor).all()
    if receptors:
        rec_branch = root_tree.add("[bold cyan]Raw Receptors & Intermediates[/bold cyan]")
        for r in receptors:
            r_color = "green" if r.status == "ready" else "yellow"
            parent_info = f" [dim](parent: {r.parent_name})[/dim]" if r.parent_name else ""
            rec_branch.add(f"[{r_color}]{r.name}[/{r_color}] [dim]({r.prep_mode})[/dim]{parent_info}")

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
