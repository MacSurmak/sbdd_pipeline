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


from cli.core.db import init_db, Receptor, Ligand, Pose, System, Trajectory


def add_receptor(args):
    """
    Registers a receptor and automatically sets up its canonical System and crystal Pose if ligand is present.
    """
    session = init_db()

    target_path = Path(args.path).resolve()
    if not target_path.exists():
        console.print(f"[bold red]FATAL: File does not exist: {target_path}[/bold red]")
        sys.exit(1)

    name = args.name
    prep_mode = args.prep
    cofactors = args.cofactors if args.cofactors else []
    ligand_name = args.ligand.strip() if args.ligand else None
    status = "ready" if prep_mode == "ready" else "unprepared"

    # 1. Register Receptor
    rec = session.query(Receptor).filter_by(name=name).first()
    if not rec:
        rec = Receptor(name=name, path=str(target_path), status=status, prep_mode=prep_mode)
        rec.meta_dict = {"cofactors": cofactors, "native_ligand": ligand_name or ""}
        session.add(rec)
        console.print(f"[bold green]✓ Registered Receptor '{name}'[/bold green] [dim]({prep_mode})[/dim]")
    else:
        rec.path = str(target_path)
        rec.status = status
        rec.prep_mode = prep_mode
        rec.meta_dict = {"cofactors": cofactors, "native_ligand": ligand_name or ""}

    # 2. Build initial System and Pose
    if ligand_name:
        lig_code = ligand_name[:4].upper()
        # Ensure Ligand entry exists
        lig = session.query(Ligand).filter_by(name=lig_code).first()
        if not lig:
            lig = Ligand(name=lig_code, path=str(target_path), status=status)
            session.add(lig)

        # Register Crystal Pose
        pose_name = f"{name}_{lig_code}_pose"
        pose = session.query(Pose).filter_by(name=pose_name).first()
        if not pose:
            pose = Pose(
                name=pose_name,
                ligand_name=lig_code,
                receptor_name=name,
                origin="crystal",
                path=str(target_path),
                status=status,
            )
            session.add(pose)

        # Register Holo System
        sys_name = f"{name}_holo"
        sys_obj = session.query(System).filter_by(name=sys_name).first()
        if not sys_obj:
            sys_obj = System(name=sys_name, receptor_name=name, pose_name=pose_name, status=status)
            sys_obj.meta_dict = {"cofactors": cofactors}
            session.add(sys_obj)
            console.print(f"[bold green]✓ Created Holo System '{sys_name}'[/bold green] [dim](Receptor: {name} + Pose: {pose_name})[/dim]")
    else:
        # Register Apo System
        sys_name = f"{name}_apo"
        sys_obj = session.query(System).filter_by(name=sys_name).first()
        if not sys_obj:
            sys_obj = System(name=sys_name, receptor_name=name, pose_name=None, status=status)
            sys_obj.meta_dict = {"cofactors": cofactors}
            session.add(sys_obj)
            console.print(f"[bold green]✓ Created Apo System '{sys_name}'[/bold green] [dim](Receptor: {name}, No Ligand)[/dim]")

    session.commit()
    session.close()


def add_ligand(args):
    """Registers unpositioned chemical molecules (single or streaming multi-molecule SDF)."""
    session = init_db()
    target_path = Path(args.path).resolve()
    if not target_path.exists():
        console.print(f"[bold red]FATAL: File does not exist: {target_path}[/bold red]")
        sys.exit(1)

    name = args.name
    existing = session.query(Ligand).filter_by(name=name).first()
    if existing:
        existing.path = str(target_path)
        existing.status = "unprepared"
        console.print(f"[yellow]Updated Ligand '{name}' in database.[/yellow]")
    else:
        new_lig = Ligand(name=name, path=str(target_path), status="unprepared")
        session.add(new_lig)
        console.print(f"[bold green]✓ Registered Ligand '{name}'[/bold green]")

    session.commit()
    session.close()


def add_traj(args):
    """Registers an MD trajectory directly linked to a physical System."""
    session = init_db()

    xtc_path = Path(args.xtc).resolve()
    tpr_path = Path(args.tpr).resolve()
    if not xtc_path.exists() or not tpr_path.exists():
        console.print("[bold red]FATAL: XTC or TPR file not found![/bold red]")
        sys.exit(1)

    target_sys = args.system
    sys_obj = session.query(System).filter_by(name=target_sys).first()
    if not sys_obj:
        # Fallback: attempt lookup by receptor name
        sys_obj = session.query(System).filter_by(receptor_name=target_sys).first()

    if not sys_obj:
        console.print(
            f"[bold red]FATAL: Target System '{target_sys}' not found in DB![/bold red]"
        )
        sys.exit(1)

    sys_name = sys_obj.name

    rep_idx = 1
    for part in args.name.split("_"):
        if part.startswith("rep") or part.startswith("R"):
            try:
                rep_idx = int("".join(filter(str.isdigit, part)))
            except ValueError:
                pass

    # Auto-detect if registering pre-cleaned trajectory
    ref_pdb = xtc_path.parent / "complex_ref.pdb"
    is_clean = xtc_path.name == "trajectory_clean.xtc" and ref_pdb.exists()

    new_traj = Trajectory(
        name=args.name,
        system_name=sys_name,
        replicate=rep_idx,
        path_xtc=str(xtc_path),
        path_tpr=str(tpr_path),
        path_clean_xtc=str(xtc_path) if is_clean else None,
        path_ref_pdb=str(ref_pdb) if is_clean else None,
        length_ns=0.0,
        status="cleaned" if is_clean else "raw",
    )
    session.merge(new_traj)
    session.commit()
    session.close()
    console.print(
        f"[bold green]✓ Registered Trajectory '{args.name}' linked to System '{sys_name}'[/bold green]"
    )
