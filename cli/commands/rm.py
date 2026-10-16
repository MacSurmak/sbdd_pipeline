"""
File: cli/commands/rm.py
Description: CLI handler to remove items from the clean 3NF database with cascade deletion.
"""

from rich.console import Console
from cli.core.db import (
    init_db,
    Target,
    ReceptorFrame,
    Substance,
    Pose,
    System,
    Trajectory,
)

console = Console()


def remove_entity(name: str) -> None:
    session = init_db()
    deleted = False

    # Check across models in hierarchical order
    target = session.query(Target).filter_by(name=name).first()
    if target:
        session.delete(target)
        console.print(
            f"[bold green]✓ Deleted Target '{name}' and all associated frames/systems (CASCADE).[/bold green]"
        )
        deleted = True

    if not deleted:
        sys_obj = session.query(System).filter_by(name=name).first()
        if sys_obj:
            session.delete(sys_obj)
            console.print(
                f"[bold green]✓ Deleted System '{name}' and all associated trajectories.[/bold green]"
            )
            deleted = True

    if not deleted:
        frame = session.query(ReceptorFrame).filter_by(name=name).first()
        if frame:
            session.delete(frame)
            console.print(
                f"[bold green]✓ Deleted ReceptorFrame '{name}' and associated poses/systems.[/bold green]"
            )
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
        sub = (
            session.query(Substance)
            .filter((Substance.name == name) | (Substance.inchikey == name))
            .first()
        )
        if sub:
            session.delete(sub)
            console.print(
                f"[bold green]✓ Deleted Chemical Substance '{name}'.[/bold green]"
            )
            deleted = True

    if deleted:
        session.commit()
    else:
        console.print(f"[bold red]✗ Item '{name}' not found in database.[/bold red]")

    session.close()
