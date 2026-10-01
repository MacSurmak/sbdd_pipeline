"""
File: cli/commands/rm.py
Description: CLI handler to remove items from the SBDD database.
"""

from rich.console import Console
from cli.core.db import init_db, Receptor, Ligand, Trajectory, PoseSet, Complex

console = Console()


def cmd_rm(args):
    session = init_db()
    target = args.target

    # Search and delete across all tables
    deleted = False

    sys_obj = session.query(System).filter_by(name=target).first()
    if sys_obj:
        session.delete(sys_obj)
        console.print(
            f"[bold green]✓ Deleted System '{target}' from database.[/bold green]"
        )
        deleted = True

    rec = session.query(Receptor).filter_by(name=target).first()
    if rec:
        session.delete(rec)
        console.print(
            f"[bold green]✓ Deleted Receptor '{target}' from database.[/bold green]"
        )
        deleted = True

    pose = session.query(Pose).filter_by(name=target).first()
    if pose:
        session.delete(pose)
        console.print(
            f"[bold green]✓ Deleted Pose '{target}' from database.[/bold green]"
        )
        deleted = True

    lig = session.query(Ligand).filter_by(name=target).first()
    if lig:
        session.delete(lig)
        console.print(
            f"[bold green]✓ Deleted Ligand '{target}' from database.[/bold green]"
        )
        deleted = True

    traj = session.query(Trajectory).filter_by(name=target).first()
    if traj:
        session.delete(traj)
        console.print(
            f"[bold green]✓ Deleted Trajectory '{target}' from database.[/bold green]"
        )
        deleted = True

    pset = session.query(PoseSet).filter_by(name=target).first()
    if pset:
        session.delete(pset)
        console.print(
            f"[bold green]✓ Deleted PoseSet '{target}' from database.[/bold green]"
        )
        deleted = True

    cmplx = session.query(Complex).filter_by(name=target).first()
    if cmplx:
        session.delete(cmplx)
        console.print(
            f"[bold green]✓ Deleted Complex '{target}' from database.[/bold green]"
        )
        deleted = True

    if deleted:
        session.commit()
    else:
        console.print(
            f"[bold red]✗ Item '{target}' not found in any database table.[/bold red]"
        )

    session.close()
