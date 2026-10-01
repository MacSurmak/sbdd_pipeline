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
        rec = Receptor(
            name=name, path=str(target_path), status=status, prep_mode=prep_mode
        )
        rec.meta_dict = {"cofactors": cofactors, "native_ligand": ligand_name or ""}
        session.add(rec)
        console.print(
            f"[bold green]✓ Registered Receptor '{name}'[/bold green] [dim]({prep_mode})[/dim]"
        )
    else:
        rec.path = str(target_path)
        rec.status = status
        rec.prep_mode = prep_mode
        rec.meta_dict = {"cofactors": cofactors, "native_ligand": ligand_name or ""}

    # 2. Build initial System and Pose with receptor-isolated chemistry
    if ligand_name:
        lig_code = ligand_name[:4].upper()
        lig_entity_name = f"{name}_{lig_code}"

        # Register Ligand uniquely bound to its structural source
        lig = session.query(Ligand).filter_by(name=lig_entity_name).first()
        if not lig:
            lig = Ligand(name=lig_entity_name, path=str(target_path), status=status)
            lig.meta_dict = {"resname": lig_code}
            session.add(lig)
        else:
            lig.path = str(target_path)
            lig.status = status

        # Register Crystal Pose
        pose_name = f"{name}_{lig_code}_pose"
        pose = session.query(Pose).filter_by(name=pose_name).first()
        if not pose:
            pose = Pose(
                name=pose_name,
                ligand_name=lig_entity_name,
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
            sys_obj = System(
                name=sys_name, receptor_name=name, pose_name=pose_name, status=status
            )
            sys_obj.meta_dict = {"cofactors": cofactors}
            session.add(sys_obj)
            console.print(
                f"[bold green]✓ Created Holo System '{sys_name}'[/bold green] [dim](Receptor: {name} + Pose: {pose_name})[/dim]"
            )
    else:
        # Register Apo System
        sys_name = f"{name}_apo"
        sys_obj = session.query(System).filter_by(name=sys_name).first()
        if not sys_obj:
            sys_obj = System(
                name=sys_name, receptor_name=name, pose_name=None, status=status
            )
            sys_obj.meta_dict = {"cofactors": cofactors}
            session.add(sys_obj)
            console.print(
                f"[bold green]✓ Created Apo System '{sys_name}'[/bold green] [dim](Receptor: {name}, No Ligand)[/dim]"
            )

    session.commit()
    session.close()


def add_ligand(args):
    """Registers unpositioned chemical molecules (single or streaming multi-molecule SDF)."""
    session = init_db()
    target_path = Path(args.path).resolve()
    if not target_path.exists():
        console.print(f"[bold red]FATAL: File does not exist: {target_path}[/bold red]")
        sys.exit(1)

    # Check molecule count by counting delimiters
    mol_count = 0
    with open(target_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("$$$$"):
                mol_count += 1

    if mol_count <= 1:
        name = args.name or target_path.stem
        smiles = None
        try:
            from rdkit import Chem

            suppl = Chem.SDMolSupplier(str(target_path), removeHs=False)
            mol = next((m for m in suppl if m is not None), None)
            if mol:
                smiles = Chem.MolToSmiles(mol)
                if (
                    not args.name
                    and mol.HasProp("_Name")
                    and mol.GetProp("_Name").strip()
                ):
                    name = mol.GetProp("_Name").strip()
        except Exception:
            pass

        existing = session.query(Ligand).filter_by(name=name).first()
        if existing:
            existing.path = str(target_path)
            existing.smiles = smiles
            existing.status = "unprepared"
            console.print(f"[yellow]Updated Ligand '{name}' in database.[/yellow]")
        else:
            new_lig = Ligand(
                name=name, path=str(target_path), smiles=smiles, status="unprepared"
            )
            session.add(new_lig)
            console.print(f"[bold green]✓ Registered Ligand '{name}'[/bold green]")
    else:
        console.print(
            f"[bold cyan]Scanning and importing library with {mol_count} molecules from {target_path.name}...[/bold cyan]"
        )
        try:
            from rdkit import Chem

            suppl = Chem.ForwardSDMolSupplier(str(target_path), removeHs=False)
            batch = []
            seen_names = set()
            idx = 1
            for mol in suppl:
                if mol is None:
                    continue
                title = (
                    mol.GetProp("_Name").strip()
                    if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
                    else f"MOL_{idx:06d}"
                )
                unique_name = title
                suffix = 1
                while (
                    unique_name in seen_names
                    or session.query(Ligand).filter_by(name=unique_name).first()
                ):
                    unique_name = f"{title}_{suffix}"
                    suffix += 1
                seen_names.add(unique_name)

                smi = Chem.MolToSmiles(mol)
                batch.append(
                    Ligand(
                        name=unique_name,
                        path=str(target_path),
                        smiles=smi,
                        status="unprepared",
                    )
                )
                idx += 1
                if len(batch) >= 1000:
                    session.bulk_save_objects(batch)
                    session.commit()
                    batch = []

            if batch:
                session.bulk_save_objects(batch)
                session.commit()

            console.print(
                f"[bold green]✓ Successfully imported {idx - 1} individual ligands into database.[/bold green]"
            )
        except ImportError:
            name = args.name or target_path.stem
            new_lig = Ligand(name=name, path=str(target_path), status="unprepared")
            session.add(new_lig)
            console.print(
                f"[bold green]✓ Registered Ligand Library '{name}' ({mol_count} molecules)[/bold green]"
            )

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

    # Strict explicit paths: only cleaned if both files are explicitly provided
    clean_xtc_path = None
    ref_pdb_path = None
    status = "raw"

    if getattr(args, "clean_xtc", None) and getattr(args, "ref_pdb", None):
        c_xtc = Path(args.clean_xtc).resolve()
        r_pdb = Path(args.ref_pdb).resolve()
        if not c_xtc.exists() or not r_pdb.exists():
            console.print(
                f"[bold red]FATAL: Explicitly passed --clean-xtc ({c_xtc}) or --ref-pdb ({r_pdb}) does not exist![/bold red]"
            )
            sys.exit(1)
        clean_xtc_path = str(c_xtc)
        ref_pdb_path = str(r_pdb)
        status = "cleaned"

    existing = session.query(Trajectory).filter_by(name=args.name).first()
    if existing:
        existing.system_name = sys_name
        existing.replicate = rep_idx
        existing.path_xtc = str(xtc_path)
        existing.path_tpr = str(tpr_path)
        existing.path_clean_xtc = clean_xtc_path
        existing.path_ref_pdb = ref_pdb_path
        existing.status = status
        console.print(
            f"[yellow]Updated Trajectory '{args.name}' linked to System '{sys_name}' ({status})[/yellow]"
        )
    else:
        traj_len = 0.0
        try:
            import MDAnalysis as mda
            u_check = mda.Universe(ref_pdb_path or str(tpr_path), clean_xtc_path or str(xtc_path))
            traj_len = round(float(u_check.trajectory.totaltime / 1000.0), 2)
        except Exception:
            pass

        new_traj = Trajectory(
            name=args.name,
            system_name=sys_name,
            replicate=rep_idx,
            path_xtc=str(xtc_path),
            path_tpr=str(tpr_path),
            path_clean_xtc=clean_xtc_path,
            path_ref_pdb=ref_pdb_path,
            length_ns=traj_len,
            status=status,
        )
        session.add(new_traj)
        console.print(
            f"[bold green]✓ Registered Trajectory '{args.name}' linked to System '{sys_name}' ({status})[/bold green]"
        )

    session.commit()
    session.close()
