"""
File: cli/commands/add.py
Description: Clean 3NF CLI handlers for importing targets, structural frames, ligands, and trajectories.
Enforces InChIKey deduplication, atomic poses, and multi-model frame transparent splitting.
"""

from pathlib import Path
import sys
import shutil
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Substance, Pose, System, Trajectory
from cli.core.chem import get_mol_identity

console = Console()


import math


def _detect_gaps(struct_path: Path) -> bool:
    """Rigorous gap detection supporting both PDB and mmCIF via Bio.PDB. Strictly forbids silent fallbacks."""
    from Bio.PDB import PDBParser, MMCIFParser, is_aa

    suffix = struct_path.suffix.lower()
    parser = (
        MMCIFParser(QUIET=True)
        if suffix in [".cif", ".mmcif"]
        else PDBParser(QUIET=True)
    )
    st = parser.get_structure("struct", str(struct_path))

    for model in st:
        for chain in model:
            prev_ca = None
            prev_res = None
            for res in chain:
                if is_aa(res) and "CA" in res:
                    res_id = res.get_id()
                    res_seq = res_id[1]
                    coord = res["CA"].get_coord()
                    if prev_ca is not None and prev_res is not None:
                        d2 = sum((coord[k] - prev_ca[k]) ** 2 for k in range(3))
                        if (res_seq - prev_res > 1) or d2 > 20.25:
                            return True
                    prev_ca = coord
                    prev_res = res_seq
    return False


STANDARD_PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "CYX",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASH",
    "GLH",
    "ACE",
    "NME",
    "NMA",
    "MSE",
    "CSE",
}


class CleanProteinSelect:
    """Strict whitelist filter: retains ONLY standard protein amino acid residues for ReceptorFrame."""

    def accept_model(self, model):
        return 1

    def accept_chain(self, chain):
        return 1

    def accept_atom(self, atom):
        return 1

    def accept_residue(self, res):
        from Bio.PDB import is_aa

        r_name = res.get_resname().strip().upper()
        if r_name in STANDARD_PROTEIN_RESIDUES or (
            is_aa(res) and res.get_id()[0] == " "
        ):
            return 1
        return 0


def _sanitize_and_split_models(
    struct_path: Path,
    dest_dir: Path,
    base_name: str,
    allowed_ligand: str | None,
    allowed_cofactors: list[str],
    allowed_ions: list[str],
) -> list[tuple[str, Path, Path]]:
    """
    Parses macromolecular structure, executes strict Fail-Fast verification of explicitly
    requested ligand, cofactors, and ions, and separates pure APO protein from raw model geometry.
    Returns list of (frame_name, apo_frame_pdb, raw_model_pdb).
    """
    from Bio.PDB import PDBParser, MMCIFParser, PDBIO

    dest_dir.mkdir(parents=True, exist_ok=True)
    suffix = struct_path.suffix.lower()

    if suffix in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True)
    elif suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"FATAL: Unsupported structure format '{suffix}'. Expected .pdb or .cif/.mmcif")

    st = parser.get_structure("struct", str(struct_path))

    # Strict Fail-Fast verification of explicitly requested molecules
    all_resnames = set()
    for model in st:
        for chain in model:
            for res in chain:
                all_resnames.add(res.get_resname().strip().upper())

    if allowed_ligand and allowed_ligand.upper() not in all_resnames:
        raise ValueError(
            f"FATAL: Explicitly specified ligand '{allowed_ligand}' was not found in {struct_path.name}! "
            f"Available residues: {sorted(list(all_resnames))}"
        )

    for cof in allowed_cofactors:
        if cof.upper() not in all_resnames:
            raise ValueError(
                f"FATAL: Explicitly specified cofactor '{cof}' was not found in {struct_path.name}! "
                f"Available residues: {sorted(list(all_resnames))}"
            )

    for ion in allowed_ions:
        if ion.upper() not in all_resnames:
            raise ValueError(
                f"FATAL: Explicitly specified ion '{ion}' was not found in {struct_path.name}! "
                f"Available residues: {sorted(list(all_resnames))}"
            )

    io = PDBIO()
    protein_selector = CleanProteinSelect()

    models = list(st.get_models())
    results = []

    if len(models) > 1:
        console.print(
            f"[cyan]Detected {len(models)} models in ensemble. Separating clean apo frames...[/cyan]"
        )
        for idx, model in enumerate(models, start=1):
            f_name = f"{base_name}_m{idx}"
            out_apo = dest_dir / f"{f_name}.pdb"
            out_raw = dest_dir / f"{f_name}_source.pdb"
            
            io.set_structure(model)
            io.save(str(out_raw))  # Preserves all atoms for heterogen extraction
            io.save(str(out_apo), protein_selector)  # Strictly APO protein
            results.append((f_name, out_apo, out_raw))
    else:
        f_name = base_name
        out_apo = dest_dir / f"{f_name}.pdb"
        out_raw = dest_dir / f"{f_name}_source.pdb"

        io.set_structure(st[0])
        io.save(str(out_raw))
        io.save(str(out_apo), protein_selector)
        results.append((f_name, out_apo, out_raw))

    return results


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


def add_receptor(
    paths: list[str],
    target: str | None = None,
    name: str | None = None,
    cofactors: list[str] | None = None,
    ions: list[str] | None = None,
    ligand: str | None = None,
) -> None:
    """
    Imports macromolecular structure(s) (PDB/CIF) into a Target.
    Enforces 100% APO protein frames and extracts explicit heterogens into atomic 3D poses.
    """
    session = init_db()
    cwd = Path.cwd()
    target_cli_name = target

    raw_inputs = paths if isinstance(paths, list) else [paths]
    struct_paths = [Path(p).resolve() for p in raw_inputs]

    for struct_path in struct_paths:
        if not struct_path.exists():
            console.print(
                f"[bold red]FATAL: Structure file not found: {struct_path}[/bold red]"
            )
            sys.exit(1)

        target_name = target_cli_name if target_cli_name else struct_path.stem
        target_obj = session.query(Target).filter_by(name=target_name).first()
        if not target_obj:
            target_obj = Target(name=target_name, description=f"Target {target_name}")
            session.add(target_obj)
            session.flush()
            console.print(f"[bold green]✓ Created Target '{target_name}'[/bold green]")

        base_name = name if (len(struct_paths) == 1 and name) else struct_path.stem
        safe_cofactors = [c.strip().upper() for c in (cofactors or [])]
        safe_ions = [i.strip().upper() for i in (ions or [])]
        lig_code = ligand.strip().upper() if ligand else None
        source_type = "crystal"

        frames_dir = cwd / "targets" / target_name / "frames"
        frames = _sanitize_and_split_models(
            struct_path, frames_dir, base_name, lig_code, safe_cofactors, safe_ions
        )

        for frame_name, apo_file, raw_file in frames:
            has_gaps = _detect_gaps(apo_file)
            initial_state = "broken" if has_gaps else "raw"

            ex_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
            if ex_frame:
                ex_frame.path = str(apo_file.resolve())
                ex_frame.state = initial_state
                frame_obj = ex_frame
                console.print(
                    f"[yellow]Updated Frame '{frame_name}' in Target '{target_name}' (state: {initial_state})[/yellow]"
                )
            else:
                frame_obj = ReceptorFrame(
                    target_id=target_obj.id,
                    name=frame_name,
                    source_type=source_type,
                    state=initial_state,
                    path=str(apo_file.resolve()),
                )
                session.add(frame_obj)
                session.flush()

                status_color = "red" if has_gaps else "green"
                console.print(
                    f"[bold green]✓ Added Frame '{frame_name}'[/bold green] [dim]([{status_color}]{initial_state}[/{status_color}], {source_type})[/dim]"
                )

            if has_gaps:
                console.print(
                    f"[!] [bold yellow]Warning: Internal missing loops detected in '{frame_name}'.[/bold yellow]\n"
                    f"    It cannot be processed by 'sbdd prep'. Provide FASTA and run 'sbdd heal' to physically repair it."
                )

            sys_name = f"{frame_name}_{lig_code}" if lig_code else frame_name
            system_poses_dir = cwd / "systems" / sys_name / "poses"
            system_poses_dir.mkdir(parents=True, exist_ok=True)

            def _process_heterogen(resname: str, sub_type: str):
                pose_sdf = system_poses_dir / f"{frame_name}_{resname}_crystal.sdf"
                _extract_ligand_atoms(raw_file, resname, pose_sdf)

                parent_inchikey, parent_smiles, mol_name, formal_charge = get_mol_identity(pose_sdf)

                sub_dir = cwd / "substances" / parent_inchikey
                sub_dir.mkdir(parents=True, exist_ok=True)
                base_sdf = sub_dir / "base.sdf"
                if not base_sdf.exists():
                    shutil.copy2(pose_sdf, base_sdf)

                sub_entity = session.query(Substance).filter_by(inchikey=parent_inchikey).first()
                if not sub_entity:
                    sub_entity = Substance(
                        inchikey=parent_inchikey,
                        smiles=parent_smiles,
                        name=mol_name or resname,
                        substance_type=sub_type,
                        state="raw",
                        batch_id="crystal_extract",
                    )
                    session.add(sub_entity)
                    session.flush()
                    console.print(
                        f"[bold green]✓ Registered {sub_type.capitalize()} '{sub_entity.name}'[/bold green] [dim](raw, {parent_inchikey})[/dim]"
                    )
                else:
                    console.print(
                        f"[bold green]✓ Matched existing {sub_type.capitalize()} '{sub_entity.name}' via InChIKey[/bold green] [dim]({parent_inchikey})[/dim]"
                    )

                pose_name = f"{frame_name}_{resname}_pose"
                p_obj = session.query(Pose).filter_by(name=pose_name).first()
                if p_obj:
                    p_obj.path = str(pose_sdf.resolve())
                    p_obj.formal_charge = formal_charge
                else:
                    p_obj = Pose(
                        substance_id=sub_entity.id,
                        receptor_frame_id=frame_obj.id,
                        name=pose_name,
                        origin="crystal",
                        state="raw",
                        formal_charge=formal_charge,
                        path=str(pose_sdf.resolve()),
                    )
                    session.add(p_obj)
                    session.flush()
                    ch_label = f"q={formal_charge:+d}" if formal_charge != 0 else "neutral"
                    console.print(
                        f"[bold green]✓ Registered Pose '{pose_name}'[/bold green] [dim]({ch_label}) -> {pose_sdf.name}[/dim]"
                    )
                return p_obj

            main_ligand_pose = None
            if lig_code:
                main_ligand_pose = _process_heterogen(lig_code, "ligand")

            hetero_poses = []
            for cof in safe_cofactors:
                c_pose = _process_heterogen(cof, "cofactor")
                hetero_poses.append(c_pose)

            for ion in safe_ions:
                i_pose = _process_heterogen(ion, "ion")
                hetero_poses.append(i_pose)

            # Purge intermediate extraction raw file to prevent disk clutter
            raw_file.unlink(missing_ok=True)

            # System creation
            ex_sys = session.query(System).filter_by(name=sys_name).first()
            if ex_sys:
                ex_sys.pose_id = main_ligand_pose.id if main_ligand_pose else None
                ex_sys.state = "raw"
                ex_sys.cofactors = hetero_poses
            else:
                new_sys = System(
                    target_id=target_obj.id,
                    receptor_frame_id=frame_obj.id,
                    pose_id=main_ligand_pose.id if main_ligand_pose else None,
                    name=sys_name,
                    state="raw",
                    origin="crystal",
                )
                new_sys.cofactors = hetero_poses
                session.add(new_sys)
                sys_type = "Holo" if main_ligand_pose else "Apo"
                console.print(
                    f"[bold green]✓ Created {sys_type} System '{sys_name}'[/bold green] [dim]with {len(hetero_poses)} cofactor(s)/ion(s)[/dim]"
                )

    session.commit()
    session.close()


def add_fasta(path: str) -> None:
    """Imports FASTA reference sequences for Targets."""
    session = init_db()
    fasta_path = Path(path).resolve()

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
                    console.print(
                        f"[yellow]Warning: Target '{current_target}' not found in DB, skipping sequence.[/yellow]"
                    )

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
    console.print(
        f"[bold green]✓ Updated reference sequences for {updated} target(s).[/bold green]"
    )


def add_ligand(path: str, name: str | None = None, force: bool = False) -> None:
    """
    Imports chemical substances (SDF/SMILES) into the Substance table.
    Performs topological deduplication via Standard InChIKey.
    """
    session = init_db()
    target_path = Path(path).resolve()
    if not target_path.exists():
        console.print(f"[bold red]FATAL: File does not exist: {target_path}[/bold red]")
        sys.exit(1)

    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:
        console.print(
            "[bold red]FATAL: RDKit is required for ligand import.[/bold red]"
        )
        sys.exit(1)

    raw_mols = [m for m in Chem.SDMolSupplier(str(target_path), removeHs=False) if m is not None]
    if not raw_mols:
        console.print(f"[bold red]FATAL: No valid molecules could be parsed from {target_path}[/bold red]")
        sys.exit(1)

    imported = 0
    skipped = 0
    uncharger = rdMolStandardize.Uncharger()
    cwd = Path.cwd()

    for idx, mol in enumerate(raw_mols, start=1):
        try:
            heavy_mol = Chem.RemoveHs(mol)
            parent_mol = uncharger.uncharge(heavy_mol)

            inchikey = Chem.MolToInchiKey(parent_mol)
            if not inchikey:
                raise ValueError("Could not compute InChIKey")
                
            smiles = Chem.MolToSmiles(parent_mol, canonical=True)
            if len(raw_mols) == 1 and name:
                mol_title = name.strip()
            elif mol.HasProp("_Name") and mol.GetProp("_Name").strip():
                mol_title = mol.GetProp("_Name").strip()
            elif name:
                mol_title = f"{name}_{idx:05d}"
            else:
                mol_title = f"MOL_{idx:05d}"

            # Canonical physical substance storage: substances/<inchikey>/base.sdf
            sub_dir = cwd / "substances" / inchikey
            sub_dir.mkdir(parents=True, exist_ok=True)
            base_sdf = sub_dir / "base.sdf"
            if not base_sdf.exists() or force:
                with Chem.SDWriter(str(base_sdf)) as w:
                    w.write(mol)

            ex = session.query(Substance).filter_by(inchikey=inchikey).first()
            if ex:
                if force:
                    ex.name = mol_title
                    ex.smiles = smiles
                    ex.path_sdf = str(base_sdf.resolve())
                else:
                    skipped += 1
                    continue
            else:
                new_sub = Substance(
                    inchikey=inchikey,
                    smiles=smiles,
                    name=mol_title,
                    substance_type="ligand",
                    state="raw",
                    batch_id=target_path.stem,
                    path_sdf=str(base_sdf.resolve()),
                )
                session.add(new_sub)

            imported += 1
            if imported % 500 == 0:
                session.commit()
        except Exception as e:
            console.print(
                f"[bold red]FATAL: Substance parsing failed at index {idx}: {e}[/bold red]"
            )
            sys.exit(1)

    session.commit()
    session.close()
    console.print(
        f"[bold green]✓ Imported {imported} unique chemical substance(s)[/bold green] [dim]({skipped} duplicates skipped)[/dim]"
    )


def add_traj(
    system: str,
    tpr: str,
    xtc: str | None = None,
    raw_dir: str | None = None,
    replicate: int = 1,
    name: str | None = None,
    length_ns: float = 0.0,
    clean_xtc: str | None = None,
    ref_pdb: str | None = None,
) -> None:
    """
    Registers an MD simulation replicate linked to a physical System.
    """
    session = init_db()

    sys_obj = session.query(System).filter_by(name=system).first()
    if not sys_obj:
        console.print(
            f"[bold red]FATAL: Target System '{system}' not found in DB![/bold red]"
        )
        sys.exit(1)

    if raw_dir:
        raw_dir_path = Path(raw_dir).resolve()
    elif xtc:
        raw_dir_path = Path(xtc).resolve().parent
    else:
        console.print(
            "[bold red]FATAL: Either raw_dir or xtc must be provided.[/bold red]"
        )
        sys.exit(1)

    tpr_path = Path(tpr).resolve()

    if not tpr_path.exists():
        console.print(f"[bold red]FATAL: TPR file not found: {tpr_path}[/bold red]")
        sys.exit(1)

    rep_idx = replicate or 1
    traj_name = name or f"{sys_obj.name}_rep{rep_idx}"

    clean_xtc_path_resolved = str(Path(clean_xtc).resolve()) if clean_xtc else None
    ref_pdb_path_resolved = str(Path(ref_pdb).resolve()) if ref_pdb else None
    status = "cleaned" if clean_xtc_path_resolved and ref_pdb_path_resolved else "raw"

    ex = session.query(Trajectory).filter_by(name=traj_name).first()
    if ex:
        ex.system_id = sys_obj.id
        ex.replicate = rep_idx
        ex.raw_dir = str(raw_dir_path)
        ex.path_tpr = str(tpr_path)
        ex.path_clean_xtc = clean_xtc_path_resolved
        ex.path_ref_pdb = ref_pdb_path_resolved
        ex.status = status
        console.print(f"[yellow]Updated Trajectory '{traj_name}' ({status})[/yellow]")
    else:
        new_traj = Trajectory(
            system_id=sys_obj.id,
            name=traj_name,
            replicate=rep_idx,
            length_ns=float(length_ns or 0.0),
            raw_dir=str(raw_dir_path),
            path_tpr=str(tpr_path),
            path_clean_xtc=clean_xtc_path_resolved,
            path_ref_pdb=ref_pdb_path_resolved,
            status=status,
        )
        session.add(new_traj)
        console.print(
            f"[bold green]✓ Registered Trajectory '{traj_name}' ({status}) linked to System '{sys_obj.name}'[/bold green]"
        )

    session.commit()
    session.close()
