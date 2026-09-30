"""
File: cli/commands/add.py
Description: CLI handlers for registering external files into the project database.
Enforces absolute paths and Fail-Fast validation.
"""

from pathlib import Path
import sys
from rich.console import Console
from cli.core.db import init_db, Receptor, Ligand

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
