"""
File: cli/commands/tree.py
Description: CLI handler to visualize the database lineage as a tree structure.
"""

from rich.console import Console
from rich.tree import Tree
from cli.core.db import init_db, Receptor, Ligand, Trajectory

console = Console()


def cmd_tree(args):
    session = init_db()

    root_tree = Tree("[bold blue]SBDD Project Database[/bold blue]")

    # 1. Receptors
    receptors = session.query(Receptor).all()
    if receptors:
        rec_branch = root_tree.add("[bold cyan]Receptors[/bold cyan]")

        # Helper to build recursive tree
        def add_children(parent_node, parent_name):
            children = [r for r in receptors if r.parent_name == parent_name]
            for child in children:
                color = "green" if child.status == "ready" else "yellow"
                node = parent_node.add(
                    f"[{color}]{child.name}[/{color}] [dim]({child.status}, {child.prep_mode}) -> {child.path}[/dim]"
                )

                # Check for Trajectories
                trajs = (
                    session.query(Trajectory).filter_by(receptor_name=child.name).all()
                )
                for t in trajs:
                    node.add(
                        f"[magenta]{t.name}[/magenta] [dim](MD: {t.length_ns} ns)[/dim]"
                    )

                add_children(node, child.name)

        # Start with root receptors (no parent)
        roots = [r for r in receptors if r.parent_name is None]
        for r in roots:
            color = "green" if r.status == "ready" else "yellow"
            node = rec_branch.add(
                f"[{color}]{r.name}[/{color}] [dim]({r.status}, {r.prep_mode}) -> {r.path}[/dim]"
            )

            # Check for Trajectories
            trajs = session.query(Trajectory).filter_by(receptor_name=r.name).all()
            for t in trajs:
                node.add(
                    f"[magenta]{t.name}[/magenta] [dim](MD: {t.length_ns} ns)[/dim]"
                )

            add_children(node, r.name)

    # 2. Ligands
    ligands = session.query(Ligand).all()
    if ligands:
        lig_branch = root_tree.add("[bold cyan]Ligands[/bold cyan]")
        for lig in ligands:
            tag = "[blue]Native[/blue]" if lig.is_native else "Library"
            lig_branch.add(f"{lig.name} [dim]({tag}) -> {lig.path}[/dim]")

    session.close()

    if not receptors and not ligands:
        console.print("[dim]Database is empty. Use 'sbdd add' to register files.[/dim]")
    else:
        console.print(root_tree)
