"""
File: cli/commands/add.py
Description: Clean 3NF CLI handlers for importing targets, structural frames, ligands, and trajectories.
Enforces InChIKey deduplication, atomic poses, and multi-model frame transparent splitting.
"""

from pathlib import Path
import sys
import shutil
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System, Trajectory
from cli.core.chem import get_mol_identity

console = Console()


import math

def _detect_gaps(pdb_path: Path) -> bool:
    """Fast mathematical detection of C-alpha chain breaks or missing residues."""
    prev_ca = None
    prev_res = None
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    try:
                        res_seq = int(line[22:26].strip())
                        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        if prev_ca is not None:
                            d2 = (x-prev_ca[0])**2 + (y-prev_ca[1])**2 + (z-prev_ca[2])**2
                            # If numbering jumps OR spatial distance between CA is > 4.5 Angstroms
                            if (res_seq - prev_res > 1) or d2 > 20.25:
                                return True
                        prev_ca = (x, y, z)
                        prev_res = res_seq
                    except ValueError:
                        continue
    except Exception:
        pass
    return False

def _count_and_split_models(
    struct_path: Path, dest_dir: Path, base_name: str
) -> list[tuple[str, Path]]:
    """
    Checks for NMR/multi-model structures in PDB/CIF.
    Splits transparently into individual files if multiple models exist.
    Returns a list of (frame_name, file_path).
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    suffix = struct_path.suffix.lower()

    if suffix not in [".pdb", ".cif", ".mmcif"]:
        out_f = dest_dir / f"{base_name}{suffix}"
        shutil.copy2(struct_path, out_f)
        return [(base_name, out_f)]

    # Detect models via simple fast scan
    models_content = []
    current_model = []
    in_model = False

    with open(struct_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("MODEL"):
                in_model = True
                current_model = [line]
            elif line.startswith("ENDMDL"):
                current_model.append(line)
                models_content.append("".join(current_model))
                current_model = []
                in_model = False
            elif in_model:
                current_model.append(line)

    if len(models_content) > 1:
        console.print(
            f"[cyan]Detected {len(models_content)} models in NMR ensemble. Splitting transparently...[/cyan]"
        )
        results = []
        for idx, m_text in enumerate(models_content, start=1):
            f_name = f"{base_name}_m{idx}"
            out_file = dest_dir / f"{f_name}.pdb"
            out_file.write_text(m_text)
            results.append((f_name, out_file))
        return results
    else:
        out_f = dest_dir / f"{base_name}{suffix}"
        shutil.copy2(struct_path, out_f)
        return [(base_name, out_f)]


def _extract_ligand_atoms(rec_path: Path, lig_resname: str, out_sdf: Path) -> None:
    """
    Extracts ligand residue from PDB, CIF, MMCIF or MAE into a standalone SDF.
    Perceives bond orders and aromaticity chemically via Schrodinger assignbondorders
    so that identical molecules produce the exact same topological InChIKey.
    """
    import os
    import subprocess
    import shutil

    lig_clean = lig_resname.strip().upper()
    out_sdf.parent.mkdir(parents=True, exist_ok=True)

    schrod_bin = shutil.which("schrodinger")
    if not schrod_bin and "SCHRODINGER" in os.environ:
        cand = Path(os.environ["SCHRODINGER"]) / "run"
        if cand.exists():
            schrod_bin = str(cand)

    if not schrod_bin:
        raise RuntimeError(
            "FATAL: Schrödinger runner ('schrodinger run') is required for robust bond order perception."
        )

    schrod_cmd = (
        [schrod_bin, "run", "python3"]
        if not schrod_bin.endswith("/run")
        else [schrod_bin, "python3"]
    )

    script_path = (
        Path(__file__).parent.parent.parent
        / "bin"
        / "utils"
        / "extract_ligand_bonds.py"
    )
    if not script_path.exists():
        raise FileNotFoundError(f"FATAL: Extraction script not found at {script_path}")

    cmd = schrod_cmd + [
        str(script_path.resolve()),
        str(rec_path.resolve()),
        lig_clean,
        str(out_sdf.resolve()),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)

    if res.returncode == 2:
        raise ValueError(
            f"FATAL: Ligand residue '{lig_clean}' not found in {rec_path.name}! {res.stderr.strip()}"
        )
    elif res.returncode != 0:
        raise RuntimeError(
            f"FATAL: Extraction failed with exit code {res.returncode}:\n{res.stderr}\n{res.stdout}"
        )

    if not out_sdf.exists() or out_sdf.stat().st_size == 0:
        raise ValueError(
            f"FATAL: Ligand extraction failed to produce valid SDF at {out_sdf}"
        )


def add_receptor(args):
    """
    Imports macromolecular structure(s) (PDB/CIF/MAE) into a Target.
    Supports single files, multiple paths, or shell globs.
    """
    session = init_db()
    cwd = Path.cwd()

    raw_inputs = args.path if isinstance(args.path, list) else [args.path]
    struct_paths = [Path(p).resolve() for p in raw_inputs]

    for struct_path in struct_paths:
        if not struct_path.exists():
            console.print(
                f"[bold red]FATAL: Structure file not found: {struct_path}[/bold red]"
            )
            sys.exit(1)

        target_name = args.target if getattr(args, "target", None) else struct_path.stem
        target = session.query(Target).filter_by(name=target_name).first()
        if not target:
            target = Target(name=target_name, description=f"Target {target_name}")
            session.add(target)
            session.flush()
            console.print(f"[bold green]✓ Created Target '{target_name}'[/bold green]")

        base_name = (
            args.name
            if (len(struct_paths) == 1 and getattr(args, "name", None))
            else struct_path.stem
        )
        cofactors = args.cofactors if getattr(args, "cofactors", None) else []
        lig_code = (
            args.ligand.strip().upper() if getattr(args, "ligand", None) else None
        )
        source_type = "crystal"

        frames_dir = cwd / "structures" / target_name
        frames = _count_and_split_models(struct_path, frames_dir, base_name)

        for frame_name, frame_file in frames:
            ex_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
            if ex_frame:
                ex_frame.path = str(frame_file.resolve())
                ex_frame.cofactors_list = cofactors
                frame_obj = ex_frame
                console.print(
                    f"[yellow]Updated Frame '{frame_name}' in Target '{target_name}'[/yellow]"
                )
            else:
                frame_obj = ReceptorFrame(
                    target_id=target.id,
                    name=frame_name,
                    source_type=source_type,
                    state="raw",
                    path=str(frame_file.resolve()),
                )
                frame_obj.cofactors_list = cofactors
                session.add(frame_obj)
                session.flush()
                console.print(
                    f"[bold green]✓ Added Frame '{frame_name}'[/bold green] [dim](raw, {source_type})[/dim]"
                )

            # Ligand & Pose Handling
            pose_obj = None
            if lig_code:
                poses_dir = cwd / "poses" / target_name
                pose_sdf = poses_dir / f"{frame_name}_{lig_code}_crystal.sdf"

                # Extract exact coordinates and protonation state
                _extract_ligand_atoms(frame_file, lig_code, pose_sdf)

                # Parent compound InChIKey for deduplication + exact formal charge for physics
                parent_inchikey, parent_smiles, mol_name, formal_charge = (
                    get_mol_identity(pose_sdf)
                )

                # Pure chemical deduplication by Parent Standard InChIKey
                lig_entity = (
                    session.query(Ligand).filter_by(inchikey=parent_inchikey).first()
                )
                if not lig_entity:
                    lig_entity = Ligand(
                        inchikey=parent_inchikey,
                        smiles=parent_smiles,
                        name=mol_name or lig_code,
                        state="raw",
                        batch_id="crystal_extract",
                    )
                    session.add(lig_entity)
                    session.flush()
                    console.print(
                        f"[bold green]✓ Registered Chemical Ligand '{lig_entity.name}'[/bold green] [dim](raw, {parent_inchikey})[/dim]"
                    )
                else:
                    console.print(
                        f"[bold green]✓ Matched existing Chemical Ligand '{lig_entity.name}' via Parent InChIKey[/bold green] [dim]({parent_inchikey})[/dim]"
                    )

                # Create Pose preserving exact physical charge state
                pose_name = f"{frame_name}_{lig_code}_pose"
                ex_pose = session.query(Pose).filter_by(name=pose_name).first()
                if ex_pose:
                    ex_pose.path = str(pose_sdf.resolve())
                    ex_pose.scores_dict = {"formal_charge": formal_charge}
                    pose_obj = ex_pose
                else:
                    pose_obj = Pose(
                        ligand_id=lig_entity.id,
                        receptor_frame_id=frame_obj.id,
                        name=pose_name,
                        origin="crystal",
                        state="raw",
                        path=str(pose_sdf.resolve()),
                    )
                    pose_obj.scores_dict = {"formal_charge": formal_charge}
                    session.add(pose_obj)
                    session.flush()
                    ch_label = (
                        f"q={formal_charge:+d}" if formal_charge != 0 else "neutral"
                    )
                    console.print(
                        f"[bold green]✓ Registered Pose '{pose_name}'[/bold green] [dim]({ch_label}) -> {pose_sdf.name}[/dim]"
                    )

            # System creation
            sys_name = f"{frame_name}_holo" if pose_obj else f"{frame_name}_apo"
            ex_sys = session.query(System).filter_by(name=sys_name).first()
            if ex_sys:
                ex_sys.pose_id = pose_obj.id if pose_obj else None
                ex_sys.state = "raw"
            else:
                new_sys = System(
                    target_id=target.id,
                    receptor_frame_id=frame_obj.id,
                    pose_id=pose_obj.id if pose_obj else None,
                    name=sys_name,
                    state="raw",
                    origin="crystal",
                )
                session.add(new_sys)
                sys_type = "Holo" if pose_obj else "Apo"
                console.print(
                    f"[bold green]✓ Created {sys_type} System '{sys_name}'[/bold green]"
                )

    session.commit()
    session.close()


def add_fasta(args):
    """Imports FASTA reference sequences for Targets."""
    session = init_db()
    fasta_path = Path(args.path).resolve()
    
    if not fasta_path.exists():
        console.print(f"[bold red]FATAL: FASTA file not found: {fasta_path}[/bold red]")
        sys.exit(1)

    updated = 0
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        current_target = None
        seq_lines = []
        
        def save_current():
            nonlocal updated
            if current_target and seq_lines:
                t = session.query(Target).filter_by(name=current_target).first()
                if t:
                    t.reference_sequence = "".join(seq_lines).replace(" ", "").upper()
                    updated += 1
                else:
                    console.print(f"[yellow]Warning: Target '{current_target}' not found in DB, skipping sequence.[/yellow]")

        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save_current()
                # Parse target name from fasta header (e.g. >EGFR_human)
                current_target = line[1:].split()[0].strip()
                seq_lines = []
            else:
                seq_lines.append(line)
        save_current()

    session.commit()
    session.close()
    console.print(f"[bold green]✓ Updated reference sequences for {updated} target(s).[/bold green]")

def add_ligand(args):
    """
    Imports chemical substances (SDF/SMILES) into the Ligand table.
    Performs topological deduplication via Standard InChIKey.
    """
    session = init_db()
    target_path = Path(args.path).resolve()
    if not target_path.exists():
        console.print(f"[bold red]FATAL: File does not exist: {target_path}[/bold red]")
        sys.exit(1)

    try:
        from rdkit import Chem
    except ImportError:
        console.print(
            "[bold red]FATAL: RDKit is required for ligand import.[/bold red]"
        )
        sys.exit(1)

    suppl = Chem.ForwardSDMolSupplier(str(target_path), removeHs=False)
    imported = 0
    skipped = 0

    for idx, mol in enumerate(suppl, start=1):
        if mol is None:
            continue
        try:
            inchikey = Chem.MolToInchiKey(mol)
            smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=True)
            name = (
                mol.GetProp("_Name").strip()
                if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
                else f"MOL_{idx:05d}"
            )

            ex = session.query(Ligand).filter_by(inchikey=inchikey).first()
            if ex:
                skipped += 1
                continue

            new_lig = Ligand(
                inchikey=inchikey,
                smiles=smiles,
                name=name,
                state="raw",
                batch_id=target_path.stem,
            )
            session.add(new_lig)
            imported += 1
            if imported % 500 == 0:
                session.commit()
        except Exception:
            continue

    session.commit()
    session.close()
    console.print(
        f"[bold green]✓ Imported {imported} unique chemical ligand(s)[/bold green] [dim]({skipped} duplicates skipped)[/dim]"
    )


def add_traj(args):
    """
    Registers an MD simulation replicate linked to a physical System.
    """
    session = init_db()

    target_sys = args.system
    sys_obj = session.query(System).filter_by(name=target_sys).first()
    if not sys_obj:
        console.print(
            f"[bold red]FATAL: Target System '{target_sys}' not found in DB![/bold red]"
        )
        sys.exit(1)

    raw_dir_path = Path(
        args.raw_dir if getattr(args, "raw_dir", None) else Path(args.xtc).parent
    ).resolve()
    tpr_path = Path(args.tpr).resolve()

    if not tpr_path.exists():
        console.print(f"[bold red]FATAL: TPR file not found: {tpr_path}[/bold red]")
        sys.exit(1)

    rep_idx = getattr(args, "replicate", 1) or 1
    name = args.name or f"{sys_obj.name}_rep{rep_idx}"

    clean_xtc_path = (
        str(Path(args.clean_xtc).resolve())
        if getattr(args, "clean_xtc", None)
        else None
    )
    ref_pdb_path = (
        str(Path(args.ref_pdb).resolve()) if getattr(args, "ref_pdb", None) else None
    )
    status = "cleaned" if clean_xtc_path and ref_pdb_path else "raw"

    ex = session.query(Trajectory).filter_by(name=name).first()
    if ex:
        ex.system_id = sys_obj.id
        ex.replicate = rep_idx
        ex.raw_dir = str(raw_dir_path)
        ex.path_tpr = str(tpr_path)
        ex.path_clean_xtc = clean_xtc_path
        ex.path_ref_pdb = ref_pdb_path
        ex.status = status
        console.print(f"[yellow]Updated Trajectory '{name}' ({status})[/yellow]")
    else:
        new_traj = Trajectory(
            system_id=sys_obj.id,
            name=name,
            replicate=rep_idx,
            length_ns=float(getattr(args, "length_ns", 0.0) or 0.0),
            raw_dir=str(raw_dir_path),
            path_tpr=str(tpr_path),
            path_clean_xtc=clean_xtc_path,
            path_ref_pdb=ref_pdb_path,
            status=status,
        )
        session.add(new_traj)
        console.print(
            f"[bold green]✓ Registered Trajectory '{name}' ({status}) linked to System '{sys_obj.name}'[/bold green]"
        )

    session.commit()
    session.close()
