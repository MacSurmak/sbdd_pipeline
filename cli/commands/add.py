"""
File: cli/commands/add.py
Description: CLI handlers for registering external files into the project database.
Enforces absolute paths and Fail-Fast validation.
"""

from pathlib import Path
import sys
from rich.console import Console
from cli.core.db import init_db, Receptor, Ligand, Trajectory

console = Console()


def add_receptor(args):
    session = init_db()

    # 1. Enforce Fail-Fast on paths
    target_path = Path(args.path).resolve()
    if not target_path.exists():
        console.print(
            f"[bold red]FATAL: Receptor file does not exist: {target_path}[/bold red]"
        )
        sys.exit(1)

    name = args.name
    prep_mode = args.prep
    cofactors = args.cofactors if args.cofactors else []
    native_ligand = args.ligand if args.ligand else ""

    # 2. Check if already exists
    existing = session.query(Receptor).filter_by(name=name).first()
    if existing:
        if args.force:
            existing.path = str(target_path)
            existing.prep_mode = prep_mode
            existing.status = "ready" if prep_mode == "ready" else "unprepared"
            existing.meta_dict = {
                "cofactors": cofactors,
                "native_ligand": native_ligand,
            }
            console.print(
                f"[yellow]Updated existing receptor '{name}' in database.[/yellow]"
            )
        else:
            console.print(
                f"[bold red]FATAL: Receptor '{name}' already exists. Use --force to overwrite.[/bold red]"
            )
            sys.exit(1)
    else:
        status = "ready" if prep_mode == "ready" else "unprepared"
        new_rec = Receptor(
            name=name, path=str(target_path), status=status, prep_mode=prep_mode
        )
        new_rec.meta_dict = {"cofactors": cofactors, "native_ligand": native_ligand}
        session.add(new_rec)
        console.print(
            f"[bold green]✓ Registered receptor '{name}'[/bold green] [dim]({prep_mode} mode)[/dim]"
        )

    session.commit()
    session.close()


def add_ligand(args):
    session = init_db()

    target_path = Path(args.path).resolve()
    if not target_path.exists():
        console.print(
            f"[bold red]FATAL: Ligand file does not exist: {target_path}[/bold red]"
        )
        sys.exit(1)

    name = args.name

    existing = session.query(Ligand).filter_by(name=name).first()
    if existing:
        if args.force:
            existing.path = str(target_path)
            existing.is_native = args.native
            existing.status = "unprepared"
            console.print(
                f"[yellow]Updated existing ligand '{name}' in database.[/yellow]"
            )
        else:
            console.print(
                f"[bold red]FATAL: Ligand '{name}' already exists. Use --force to overwrite.[/bold red]"
            )
            sys.exit(1)
    else:
        new_lig = Ligand(
            name=name, path=str(target_path), is_native=args.native, status="unprepared"
        )
        session.add(new_lig)
        native_tag = "[bold cyan](Native)[/bold cyan]" if args.native else "(Library)"
        console.print(
            f"[bold green]✓ Registered ligand '{name}' {native_tag}[/bold green]"
        )

    session.commit()
    session.close()


def add_traj(args):
    session = init_db()

    xtc_path = Path(args.xtc).resolve()
    tpr_path = Path(args.tpr).resolve()

    if not xtc_path.exists():
        console.print(
            f"[bold red]FATAL: XTC file does not exist: {xtc_path}[/bold red]"
        )
        sys.exit(1)
    if not tpr_path.exists():
        console.print(
            f"[bold red]FATAL: TPR file does not exist: {tpr_path}[/bold red]"
        )
        sys.exit(1)

    rec = session.query(Receptor).filter_by(name=args.receptor).first()
    if not rec:
        console.print(
            f"[bold red]FATAL: Parent Receptor '{args.receptor}' not found in DB![/bold red]"
        )
        sys.exit(1)

    name = args.name
    existing = session.query(Trajectory).filter_by(name=name).first()
    if existing:
        if args.force:
            existing.path_xtc = str(xtc_path)
            existing.path_tpr = str(tpr_path)
            existing.receptor_name = args.receptor
            console.print(
                f"[yellow]Updated existing trajectory '{name}' in database.[/yellow]"
            )
        else:
            console.print(
                f"[bold red]FATAL: Trajectory '{name}' already exists. Use --force to overwrite.[/bold red]"
            )
            sys.exit(1)
    else:
        new_traj = Trajectory(
            name=name,
            receptor_name=args.receptor,
            path_xtc=str(xtc_path),
            path_tpr=str(tpr_path),
            length_ns=0.0,
        )
        session.add(new_traj)
        console.print(
            f"[bold green]✓ Registered MD trajectory '{name}' linked to receptor '{args.receptor}'[/bold green]"
        )

    session.commit()
    session.close()
