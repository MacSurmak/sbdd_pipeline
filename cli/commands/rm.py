"""
File: cli/commands/rm.py
Description: CLI handler to remove items from the clean 3NF database with cascade deletion.
"""

from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System, Trajectory

console = Console()


def cmd_rm(args):
    session = init_db()
    name = args.target
    deleted = False

    # Check across models in hierarchical order
    target = session.query(Target).filter_by(name=name).first()
    if target:
        session.delete(target)
        console.print(f"[bold green]✓ Deleted Target '{name}' and all associated frames/systems (CASCADE).[/bold green]")
        deleted = True

    if not deleted:
        sys_obj = session.query(System).filter_by(name=name).first()
        if sys_obj:
            session.delete(sys_obj)
            console.print(f"[bold green]✓ Deleted System '{name}' and all associated trajectories.[/bold green]")
            deleted = True

    if not deleted:
        frame = session.query(ReceptorFrame).filter_by(name=name).first()
        if frame:
            session.delete(frame)
            console.print(f"[bold green]✓ Deleted ReceptorFrame '{name}' and associated poses/systems.[/bold green]")
            deleted = True

    if not deleted:
        pose = session.query(Pose).filter_by(name=name).first()
        if pose:
            session.delete(pose)
            console.print(f"[bold green]✓ Deleted Pose '{name}'.[/bold green]")
            deleted = True

    if not deleted:
        traj = session.query(Trajectory).filter_by(name=name).first()
        if traj:
            session.delete(traj)
            console.print(f"[bold green]✓ Deleted Trajectory '{name}'.[/bold green]")
            deleted = True

    if not deleted:
        lig = session.query(Ligand).filter((Ligand.name == name) | (Ligand.inchikey == name)).first()
        if lig:
            session.delete(lig)
            console.print(f"[bold green]✓ Deleted Chemical Ligand '{name}'.[/bold green]")
            deleted = True

    if deleted:
        session.commit()
    else:
        console.print(f"[bold red]✗ Item '{name}' not found in database.[/bold red]")

    session.close()
