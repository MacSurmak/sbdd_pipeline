"""
File: cli/commands/dock.py
Description: CLI handler for 'sbdd dock'. Reads DB to construct Cross-Docking plan.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import (
    init_db,
    Target,
    ReceptorFrame,
    Substance,
    SubstanceState,
    Pose,
    System,
)
from cli.core.chem import get_microstate_identity
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def _format_metric(val):
    if val is None or val == "" or str(val).lower() == "nan":
        return "-"
    try:
        return f"{float(val):.2f}"
    except (ValueError, TypeError):
        return str(val)


def _color_rmsd(val_str):
    if val_str is None or val_str == "" or str(val_str).lower() in ["nan", "-"]:
        return "-"
    try:
        v = float(val_str)
        return (
            f"[bold green]{v:.2f}[/bold green]"
            if v <= 2.0
            else f"[bold red]{v:.2f}[/bold red]"
        )
    except (ValueError, TypeError):
        return str(val_str)


def run_dock(
    target: str | None = None,
    frame_names: list[str] | None = None,
    ligand_names: list[str] | None = None,
    self_dock: bool = False,
    name: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """Run distributed docking into ReceptorFrames."""
    cwd = Path.cwd()
    session = init_db()
    is_self = self_dock

    # 1. Fetch requested Ready Frames
    query = session.query(ReceptorFrame).filter_by(state="ready")
    if target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif frame_names:
        query = query.filter(ReceptorFrame.name.in_(frame_names))

    frames = query.all()
    if not frames:
        console.print(
            "[bold red]FATAL: No ready ReceptorFrames found to dock. Run 'sbdd prep' first![/bold red]"
        )
        sys.exit(1)

    # 2. Build Plan (Self 1-to-1 vs Matrix Screening)
    plan = {
        "mode": "self" if is_self else "matrix",
        "receptors": [],
        "ligands": [],
        "pairs": [],
    }

    if is_self:
        active_frames = []

        for r in frames:
            target_sys = (
                session.query(System)
                .filter(System.receptor_frame_id == r.id, System.pose_id.isnot(None))
                .first()
            )
            if not target_sys or not target_sys.pose:
                continue

            native_pose = target_sys.pose
            native_lig = native_pose.substance

            active_frames.append((r, native_pose, native_lig))

        if not active_frames:
            console.print(
                f"[bold green]✓ All {len(frames)} ready frame(s) already have redocking poses completed. Nothing to do![/bold green]"
            )
            session.close()
            return

        for r, native_pose, native_lig in active_frames:
            # Pure 3NF Asset Resolution: zero heuristics, zero with_suffix()
            if not r.path_mae or not r.path_pdb:
                raise ValueError(
                    f"FATAL: ReceptorFrame '{r.name}' does not have both path_mae and path_pdb "
                    f"registered in the 3NF database! Entity must be prepared cleanly."
                )

            path_mae = Path(r.path_mae).resolve()
            path_pdb = Path(r.path_pdb).resolve()

            if not path_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered receptor MAE does not exist: {path_mae}"
                )
            if not path_pdb.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered receptor PDB does not exist: {path_pdb}"
                )

            if not native_pose.path_sdf or not native_pose.path_mae:
                raise ValueError(
                    f"FATAL: Native Pose for '{r.name}' does not have both path_sdf and path_mae "
                    f"registered in the 3NF database!"
                )

            pose_sdf = Path(native_pose.path_sdf).resolve()
            pose_mae = Path(native_pose.path_mae).resolve()

            if not pose_sdf.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered pose SDF does not exist: {pose_sdf}"
                )
            if not pose_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Registered pose MAE does not exist: {pose_mae}"
                )

            ion_sdfs = []
            ion_maes = []
            ion_resnames = []
            cof_sdfs = []
            cof_maes = []
            cof_resnames = []

            for c in target_sys.cofactors:
                if c.substance.substance_type == "ion":
                    ion_sdfs.append(str(Path(c.path_sdf).resolve()))
                    if c.path_mae:
                        ion_maes.append(str(Path(c.path_mae).resolve()))
                    ion_resnames.append(c.res_name)
                else:
                    cof_sdfs.append(str(Path(c.path_sdf).resolve()))
                    if c.path_mae:
                        cof_maes.append(str(Path(c.path_mae).resolve()))
                    cof_resnames.append(c.res_name)

            plan["receptors"].append(
                {
                    "name": r.name,
                    "path_mae": str(path_mae),
                    "path_pdb": str(path_pdb),
                    "ref_sdf": str(pose_sdf),
                    "ref_mae": str(pose_mae),
                    "native_ligand": native_pose.res_name or native_lig.name,
                    "ion_sdfs": ion_sdfs,
                    "ion_maes": ion_maes,
                    "ion_resnames": ion_resnames,
                    "cofactor_sdfs": cof_sdfs,
                    "cofactor_maes": cof_maes,
                    "cofactor_resnames": cof_resnames,
                }
            )

            unique_lig_name = f"{r.name}_{native_lig.name}"

            # Strict validation: Candidate ligand MUST be in state 'ready' (processed via LigPrep)
            if (
                native_lig.state != "ready"
                or not native_lig.path_mae
                or not native_lig.path_sdf
            ):
                console.print(
                    f"[bold red]FATAL: Chemical Ligand '{native_lig.name}' is in state '{native_lig.state}'! "
                    f"Run 'sbdd prep' first so LigPrep generates canonical 3D states.[/bold red]"
                )
                session.close()
                sys.exit(1)

            lig_mae = Path(native_lig.path_mae).resolve()
            lig_sdf = Path(native_lig.path_sdf).resolve()

            if not lig_mae.exists() or not lig_sdf.exists():
                console.print(
                    f"[bold red]FATAL: LigPrep deliverables missing for '{native_lig.name}': {lig_mae}. Run 'sbdd prep'.[/bold red]"
                )
                session.close()
                sys.exit(1)

            plan["ligands"].append(
                {
                    "name": unique_lig_name,
                    "path_sdf": str(lig_sdf),
                    "path_mae": str(lig_mae),
                    "is_native": True,
                    "receptor": r.name,
                    "ligand_title": native_lig.name,
                    "substance_id": native_lig.id,
                }
            )

        console.print(
            f"[bold green]✓ 1-to-1 Redocking Mode: Matched {len(active_frames)} receptors directly to their native ligands.[/bold green]"
        )

    else:
        if not ligand_names:
            console.print(
                "[bold red]FATAL: Parameter --ligands is required for cross-docking matrix mode.[/bold red]"
            )
            sys.exit(1)

        ligands = (
            session.query(Substance)
            .filter(
                Substance.name.in_(ligand_names), Substance.substance_type == "ligand"
            )
            .all()
        )
        if len(ligands) != len(ligand_names):
            found = {l.name for l in ligands}
            missing = set(ligand_names) - found
            console.print(
                f"[bold red]FATAL: Ligands not found in DB: {missing}[/bold red]"
            )
            sys.exit(1)

        for r in frames:
            sys_obj = session.query(System).filter_by(receptor_frame_id=r.id).first()
            if not sys_obj or not sys_obj.pose:
                console.print(
                    f"[bold red]FATAL: Receptor '{r.name}' has no reference pose! Screening requires a pocket reference.[/bold red]"
                )
                sys.exit(1)

            pose_obj = sys_obj.pose
            native_lig = pose_obj.substance.name

            if not r.path_mae or not r.path_pdb:
                raise ValueError(
                    f"FATAL: ReceptorFrame '{r.name}' does not have both path_mae and path_pdb registered!"
                )

            path_mae = Path(r.path_mae).resolve()
            path_pdb = Path(r.path_pdb).resolve()

            if not path_mae.exists() or not path_pdb.exists():
                raise FileNotFoundError(
                    f"FATAL: Receptor deliverables missing for '{r.name}' on disk"
                )

            if not pose_obj.path_sdf or not pose_obj.path_mae:
                raise ValueError(
                    f"FATAL: Reference Pose '{pose_obj.name}' missing path_sdf or path_mae registration!"
                )

            ref_sdf = Path(pose_obj.path_sdf).resolve()
            ref_mae = Path(pose_obj.path_mae).resolve()

            if not ref_sdf.exists() or not ref_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Reference pose deliverables missing for '{r.name}' on disk"
                )

            ion_sdfs = []
            ion_maes = []
            ion_resnames = []
            cof_sdfs = []
            cof_maes = []
            cof_resnames = []

            for c in sys_obj.cofactors:
                if c.substance.substance_type == "ion":
                    ion_sdfs.append(str(Path(c.path_sdf).resolve()))
                    if c.path_mae:
                        ion_maes.append(str(Path(c.path_mae).resolve()))
                    ion_resnames.append(c.res_name)
                else:
                    cof_sdfs.append(str(Path(c.path_sdf).resolve()))
                    if c.path_mae:
                        cof_maes.append(str(Path(c.path_mae).resolve()))
                    cof_resnames.append(c.res_name)

            plan["receptors"].append(
                {
                    "name": r.name,
                    "path_mae": str(path_mae),
                    "path_pdb": str(path_pdb),
                    "ref_sdf": str(ref_sdf),
                    "ref_mae": str(ref_mae),
                    "native_ligand": native_lig,
                    "ion_sdfs": ion_sdfs,
                    "ion_maes": ion_maes,
                    "ion_resnames": ion_resnames,
                    "cofactor_sdfs": cof_sdfs,
                    "cofactor_maes": cof_maes,
                    "cofactor_resnames": cof_resnames,
                }
            )

        for l in ligands:
            if l.state != "ready" or not l.path_mae or not l.path_sdf:
                console.print(
                    f"[bold red]FATAL: Chemical Ligand '{l.name}' is in state '{l.state}'! Run 'sbdd prep' first.[/bold red]"
                )
                session.close()
                sys.exit(1)

            lig_mae = Path(l.path_mae).resolve()
            lig_sdf = Path(l.path_sdf).resolve()

            if not lig_mae.exists() or not lig_sdf.exists():
                console.print(
                    f"[bold red]FATAL: Prepared files missing for '{l.name}' ({lig_mae}). Run 'sbdd prep' first.[/bold red]"
                )
                session.close()
                sys.exit(1)

            has_crystal_pose = (
                session.query(Pose)
                .filter_by(substance_id=l.id, origin="crystal")
                .first()
                is not None
            )

            plan["ligands"].append(
                {
                    "name": l.name,
                    "path_sdf": str(lig_sdf),
                    "path_mae": str(lig_mae),
                    "is_native": has_crystal_pose,
                    "ligand_title": l.name,
                    "substance_id": l.id,
                }
            )

        console.print(
            f"[bold cyan]Cross-Docking Matrix Plan generated: {len(frames)} Receptors × {len(ligands)} Ligand Sets[/bold cyan]"
        )

    plan_file = cwd / "dock_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan, f, indent=2)

    # 4. Build Nextflow args
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="screen")
    nf_args.extend(["--dock_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    # 5. Run execution
    ret_code = execute_nextflow(
        "screen", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow docking failed.[/bold red]")
        sys.exit(ret_code)

    # 6. Register results in Database strictly via Plan & Receipt protocol
    num_ligs = len(plan.get("ligands", []))
    default_name = (
        f"redock_{len(frames)}" if is_self else f"screen_{len(frames)}x{num_ligs}"
    )
    poses_name = name if name else default_name

    receipt_file = cwd / "dock_receipt.tsv"
    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not generated by Nextflow. "
            "Docking failed or terminated without producing deliverables."
        )

    res_csv = cwd / "reports" / "ranking.csv"
    if not res_csv.exists() or res_csv.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Ranking summary CSV {res_csv} does not exist. Pipeline failed to summarize docking."
        )

    from rdkit import Chem
    import csv
    from rich.table import Table

    # Build deterministic substance lookup from plan (zero string splitting heuristics)
    # Crucially, map by the true chemical ligand_title, not the Nextflow routing name
    substance_by_title = {
        item.get("ligand_title", item["name"]): session.query(Substance).get(
            item["substance_id"]
        )
        for item in plan["ligands"]
    }

    # Parse receipt: map of (frame_name, engine) -> absolute validated merged SDF path
    receipt_map = {}
    with open(receipt_file, "r", encoding="utf-8") as rf:
        for line_num, line in enumerate(rf, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                session.close()
                raise ValueError(
                    f"FATAL: Corrupted line {line_num} in {receipt_file}: '{line}'. Expected 3 tab-separated columns."
                )
            r_frame_id, r_engine, r_sdf_path = parts[0], parts[1], parts[2]
            sdf_file = Path(r_sdf_path).resolve()
            if not sdf_file.exists() or sdf_file.stat().st_size == 0:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: Published pose file specified in receipt does not exist or is empty: {sdf_file}"
                )
            receipt_map[(r_frame_id, r_engine)] = sdf_file

    def _sanitize_sdf_3d(path: Path):
        """Physically rewrite MDL V2000 headers to enforce 3D tags, eradicating RDKit C++ warnings on stdout."""
        try:
            with open(path, "r", encoding="utf-8") as f:
                lines = f.read().splitlines()
            changed = False
            for i in range(len(lines)):
                # Dimension tag is strictly on line 2 (index 1) of every MOL block
                if i == 1 or (i > 1 and lines[i - 2] == "$$$$"):
                    if len(lines[i]) >= 22 and lines[i][20:22] != "3D":
                        lines[i] = lines[i][:20] + "3D" + lines[i][22:]
                        changed = True
                    elif len(lines[i]) < 22:
                        lines[i] = lines[i].ljust(20) + "3D"
                        changed = True
            if changed:
                with open(path, "w", encoding="utf-8") as f:
                    f.write("\n".join(lines) + "\n")
        except Exception:
            pass

    # Index molecules in merged docking SDFs by title to slice individual pose SDFs
    poses_by_engine_and_title = {}
    for (f_name, eng), m_sdf in receipt_map.items():
        _sanitize_sdf_3d(m_sdf)
        suppl = Chem.SDMolSupplier(str(m_sdf), removeHs=False)
        for mol_idx, m in enumerate(suppl):
            if m is None:
                continue
            m_title = m.GetProp("_Name").strip() if m.HasProp("_Name") else ""
            if not m_title and m.HasProp("s_m_title"):
                m_title = m.GetProp("s_m_title").strip()
            poses_by_engine_and_title.setdefault((f_name, eng, m_title), []).append(m)

    table = Table(title="🎯 Redocking & Scoring Summary", show_header=True)
    table.add_column("Receptor", style="cyan")
    table.add_column("Ligand", style="bold")
    table.add_column("GlideScore", justify="right")
    table.add_column("CNN Affinity", justify="right")
    table.add_column("MM-GBSA (kcal/mol)", justify="right")
    table.add_column("Glide RMSD (Å)", justify="center")
    table.add_column("GNINA RMSD (Å)", justify="center")
    table.add_column("Verdict", justify="center")

    poses_count = 0
    pose_counter = {}

    with open(res_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("Title", "").strip()
            frame_name = row.get("Receptor", row.get("Frame", "")).strip()
            if not title or not frame_name:
                session.close()
                raise ValueError(
                    f"FATAL: Corrupted ranking row with missing Title or Frame: {row}"
                )

            r_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
            if not r_frame:
                session.close()
                raise ValueError(
                    f"FATAL: ReceptorFrame '{frame_name}' in ranking is not registered in SQLite database!"
                )

            lig = substance_by_title.get(title)
            if not lig:
                session.close()
                raise KeyError(
                    f"FATAL: Ligand Title '{title}' in ranking was not declared in docking plan!"
                )

            engines_for_frame = [
                eng for (f_id, eng) in receipt_map.keys() if f_id == frame_name
            ]
            if not engines_for_frame:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: No receipt records found for receptor frame '{frame_name}' in {receipt_file}!"
                )

            for eng in engines_for_frame:
                mol_pool = poses_by_engine_and_title.get((frame_name, eng, title), [])
                if not mol_pool:
                    session.close()
                    raise ValueError(
                        f"FATAL: 3D coordinates for molecule '{title}' not found in {eng} docking deliverable!"
                    )

                # Consume next 3D conformation for this ligand/engine pair
                mol_conformation = mol_pool.pop(0)

                # Topological Restoration for GNINA/United-Atom engines
                # We strictly graft the new Heavy Atom coordinates onto the perfect LigPrep template
                # to recover lost hydrogens, bond orders, and formal charges before Antechamber processing.
                if eng != "glide":
                    try:
                        from rdkit.Chem import AllChem
                        # Load the canonical SSOT template prepared by LigPrep
                        ref_mol = next(Chem.SDMolSupplier(str(Path(lig.path_sdf).resolve()), removeHs=False), None)
                        if ref_mol:
                            # 1. Strip hydrogens from both to establish a 1:1 heavy-atom mapping isomorphism
                            ref_heavy = Chem.RemoveHs(ref_mol)
                            dock_heavy = Chem.RemoveHs(mol_conformation)
                            match = ref_heavy.GetSubstructMatch(dock_heavy)

                            if match and len(match) == dock_heavy.GetNumAtoms():
                                # 2. Graft coordinates onto the perfect heavy-atom template
                                conf = ref_heavy.GetConformer()
                                dock_conf = dock_heavy.GetConformer()
                                for dock_idx, ref_idx in enumerate(match):
                                    conf.SetAtomPosition(ref_idx, dock_conf.GetAtomPosition(dock_idx))

                                # 3. Reconstruct ideal hydrogen coordinates based on the new heavy atom geometry
                                restored_mol = Chem.AddHs(ref_heavy, addCoords=True)
                                # Transfer properties and title
                                restored_mol.SetProp("_Name", mol_conformation.GetProp("_Name") if mol_conformation.HasProp("_Name") else "")
                                mol_conformation = restored_mol
                    except Exception as e:
                        console.print(f"[yellow]Warning: Failed to topologically restore GNINA pose '{title}': {e}. Storing as-is.[/yellow]")

                key = (frame_name, title, eng)
                pose_counter[key] = pose_counter.get(key, 0) + 1
                p_idx = pose_counter[key]

                p_suffix = (
                    f"redock_{eng}_p{p_idx}" if is_self else f"dock_{eng}_p{p_idx}"
                )
                pose_name = f"{frame_name}_{title}_{p_suffix}"
                sys_name = f"{frame_name}_{title}_{p_suffix}"
                origin_tag = f"redock_{eng}" if is_self else f"docking_{eng}"

                # Canonical 3NF Domain Placement: substances/<inchikey>/poses/<pose_name>/<res>_pose.sdf
                sub_pose_dir = cwd / "substances" / lig.inchikey / "poses" / pose_name
                sub_pose_dir.mkdir(parents=True, exist_ok=True)
                clean_res = lig.res_name or title[:4].strip().upper()
                pose_sdf_file = sub_pose_dir / f"{clean_res}_pose.sdf"

                with Chem.SDWriter(str(pose_sdf_file)) as writer:
                    writer.write(mol_conformation)

                # Identify the exact chemical microstate (protomer, tautomer, stereoisomer)
                (
                    _,
                    _,
                    _,
                    state_hash,
                    isomeric_smiles,
                    formal_charge,
                ) = get_microstate_identity(mol_conformation)

                state_entity = (
                    session.query(SubstanceState)
                    .filter_by(substance_id=lig.id, state_hash=state_hash)
                    .first()
                )
                if not state_entity:
                    state_dir = cwd / "substances" / lig.inchikey / "states" / state_hash
                    state_dir.mkdir(parents=True, exist_ok=True)
                    state_sdf_file = state_dir / f"{clean_res}_state.sdf"
                    with Chem.SDWriter(str(state_sdf_file)) as writer:
                        writer.write(mol_conformation)

                    state_entity = SubstanceState(
                        substance_id=lig.id,
                        state_hash=state_hash,
                        isomeric_smiles=isomeric_smiles,
                        formal_charge=formal_charge,
                        state="raw",
                        path_sdf=str(state_sdf_file.resolve()),
                    )
                    session.add(state_entity)
                    session.flush()

                f_charge = formal_charge
                scores_dict = {}
                if eng == "glide":
                    if row.get("GlideScore") and row["GlideScore"] != "nan":
                        scores_dict["GlideScore"] = float(row["GlideScore"])
                    if (
                        row.get("InSitu_RMSD_Glide_A")
                        and row["InSitu_RMSD_Glide_A"] != "nan"
                    ):
                        scores_dict["RMSD_Glide_A"] = float(row["InSitu_RMSD_Glide_A"])
                elif eng == "gnina":
                    if row.get("Dock_CNNaffinity") and row["Dock_CNNaffinity"] != "nan":
                        scores_dict["CNNaffinity"] = float(row["Dock_CNNaffinity"])
                    if row.get("Dock_CNNscore") and row["Dock_CNNscore"] != "nan":
                        scores_dict["CNNscore"] = float(row["Dock_CNNscore"])
                    if row.get("Dock_VinaAff") and row["Dock_VinaAff"] != "nan":
                        scores_dict["VinaAff"] = float(row["Dock_VinaAff"])
                    if (
                        row.get("InSitu_RMSD_GNINA_A")
                        and row["InSitu_RMSD_GNINA_A"] != "nan"
                    ):
                        scores_dict["RMSD_GNINA_A"] = float(row["InSitu_RMSD_GNINA_A"])

                if row.get("Ultimate_Score") and row["Ultimate_Score"] != "nan":
                    scores_dict["Ultimate_Score"] = float(row["Ultimate_Score"])

                ex_pose = session.query(Pose).filter_by(name=pose_name).first()
                if ex_pose:
                    ex_pose.substance_state_id = state_entity.id
                    ex_pose.path = str(pose_sdf_file.resolve())
                    ex_pose.path_sdf = str(pose_sdf_file.resolve())
                    ex_pose.res_name = clean_res
                    ex_pose.scores_dict = scores_dict
                    ex_pose.origin = origin_tag
                    ex_pose.state = "ready"
                    pose_obj = ex_pose
                else:
                    new_p = Pose(
                        substance_state_id=state_entity.id,
                        receptor_frame_id=r_frame.id,
                        name=pose_name,
                        res_name=clean_res,
                        origin=origin_tag,
                        state="ready",
                        path=str(pose_sdf_file.resolve()),
                        path_sdf=str(pose_sdf_file.resolve()),
                    )
                    new_p.scores_dict = scores_dict
                    session.add(new_p)
                    session.flush()
                    pose_obj = new_p

                # Deterministic 3NF inheritance: query cofactors registered on this frame or lineage parent
                lookup_frame_id = (
                    r_frame.parent_frame_id if r_frame.parent_frame_id else r_frame.id
                )
                parent_cofactors = (
                    session.query(Pose)
                    .join(SubstanceState)
                    .join(Substance)
                    .filter(
                        Pose.receptor_frame_id == lookup_frame_id,
                        Pose.origin == "crystal",
                        Substance.substance_type.in_(["cofactor", "ion"]),
                    )
                    .all()
                )

                ex_sys = session.query(System).filter_by(name=sys_name).first()
                if ex_sys:
                    ex_sys.pose_id = pose_obj.id
                    ex_sys.state = "raw"
                    ex_sys.cofactors = parent_cofactors
                else:
                    new_sys = System(
                        target_id=r_frame.target_id,
                        receptor_frame_id=r_frame.id,
                        pose_id=pose_obj.id,
                        name=sys_name,
                        state="raw",
                        origin=origin_tag,
                    )
                    new_sys.cofactors = parent_cofactors
                    session.add(new_sys)

                poses_count += 1

            # Populate display row for the Rich Summary Table
            gscore = row.get("GlideScore", "-")
            cnn = row.get("Resc_CNNaffinity") or row.get("Dock_CNNaffinity") or "-"
            mmgbsa = row.get("MMGBSA", "-")
            g_rmsd = row.get("InSitu_RMSD_Glide_A", "-")
            n_rmsd = row.get("InSitu_RMSD_GNINA_A", "-")

            disp_glide = _color_rmsd(g_rmsd)
            disp_gnina = _color_rmsd(n_rmsd)

            vg = float(g_rmsd) if g_rmsd and g_rmsd not in ["nan", "-"] else None
            vn = float(n_rmsd) if n_rmsd and n_rmsd not in ["nan", "-"] else None
            if (vg is not None and vg <= 2.5) or (vn is not None and vn <= 2.5):
                verdict = "[bold green]✓ SUCCESS[/bold green]"
            elif vg is not None or vn is not None:
                verdict = "[bold yellow]⚠ DRIFT[/bold yellow]"
            else:
                verdict = "[dim]-[/dim]"

            table.add_row(
                frame_name,
                title,
                _format_metric(gscore),
                _format_metric(cnn),
                _format_metric(mmgbsa),
                disp_glide,
                disp_gnina,
                verdict,
            )

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Docking verified via receipt! Registered {poses_count} pose(s) in DB as '{poses_name}'.[/bold green]\n"
    )
    console.print(table)
