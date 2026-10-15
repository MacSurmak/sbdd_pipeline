"""
File: cli/commands/prep.py
Description: CLI handler for 'sbdd prep' using clean 3NF database lineage.
Prepares raw ReceptorFrames and in-situ Poses into ready states without mutating raw records.
"""

import json
from pathlib import Path
import sys
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Substance, Pose, System
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def run_prep(
    target: str | None = None,
    frames: list[str] | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """Orchestrates Tri-partite preparation: Receptors (PrepWizard), Substances (LigPrep), and Poses (Antechamber)."""
    cwd = Path.cwd()
    session = init_db()

    # 1. Fetch Receptors needing PrepWizard
    query = session.query(ReceptorFrame).filter(ReceptorFrame.state == "raw")
    if target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found in database.[/bold red]")
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif frames:
        query = query.filter(ReceptorFrame.name.in_(frames))

    raw_frames = query.all()

    # 2. Fetch abstract library Substances needing LigPrep (only type = 'ligand')
    query_substances = session.query(Substance).filter_by(state="raw", substance_type="ligand").all()

    # 3. Fetch in-situ Poses lacking AM1-BCC quantum parameters (excluding monoatomic ions)
    query_poses = session.query(Pose).join(Substance).filter(
        Pose.path_mol2.is_(None),
        Substance.substance_type != "ion"
    ).all()

    if not raw_frames and not query_substances and not query_poses:
        console.print("[green]No raw frames, substances, or unparametrized poses found. Everything is up to date![/green]")
        session.close()
        return

    # Build Plan
    rec_plan = []
    frame_map = {}
    for frame in raw_frames:
        frame_map[frame.name] = frame
        system = session.query(System).filter_by(receptor_frame_id=frame.id).first()
        native_lig_code = system.pose.substance.name if (system and system.pose and system.pose.substance) else ""
        cofactor_names = [c.substance.name for c in system.cofactors if c.substance] if system else []
        rec_plan.append({
            "name": frame.name,
            "path": str(Path(frame.path).resolve()),
            "prep_mode": "full",
            "cofactors": cofactor_names,
            "native_ligand": native_lig_code,
        })

    sub_plan = []
    if query_substances:
        raw_lig_dir = (cwd / "prepared" / "ligands" / "raw").resolve()
        raw_lig_dir.mkdir(parents=True, exist_ok=True)
        for sub in query_substances:
            sub_raw_sdf = raw_lig_dir / f"{sub.name}.sdf"
            
            from rdkit import Chem
            from rdkit.Chem import AllChem
            m = Chem.MolFromSmiles(sub.smiles)
            if not m:
                raise ValueError(f"FATAL: RDKit failed to parse SMILES for '{sub.name}'")
            m = Chem.AddHs(m)
            AllChem.Compute2DCoords(m)
            m.SetProp("_Name", sub.name)
            with Chem.SDWriter(str(sub_raw_sdf)) as w:
                w.write(m)

            sub_plan.append({
                "name": sub.name,
                "path": str(sub_raw_sdf),
                "status": "unprepared"
            })

    pose_plan = []
    for pose in query_poses:
        ch = pose.scores_dict.get("formal_charge", 0)
        res_name = pose.substance.name[:4].upper() if pose.substance else "LIG"
        pose_plan.append({
            "name": pose.name,
            "path": pose.path,
            "charge": ch,
            "res_name": res_name
        })

    plan_data = {"receptors": rec_plan, "ligands": sub_plan, "poses": pose_plan}
    plan_file = cwd / "prep_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_data, f, indent=2)

    console.print(f"[bold cyan]Generated prep_plan.json ({len(rec_plan)} receptors, {len(sub_plan)} substances, {len(pose_plan)} poses).[/bold cyan]")

    # Execute Nextflow
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="prepare")
    nf_args.extend(["--plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    ret_code = execute_nextflow("prepare", nf_args, cwd, profile="local" if local else "slurm", resume=resume, configs=configs)
    if ret_code != 0:
        console.print("[bold red]Nextflow preparation failed. Database not updated.[/bold red]")
        sys.exit(ret_code)

    # Read Receipts
    prep_rec_count, prep_sub_count, prep_pose_count = 0, 0, 0
    rec_receipt = cwd / "prep_receipt.tsv"
    if rec_receipt.exists():
        for line in rec_receipt.read_text().strip().splitlines():
            if not line.strip(): continue
            orig_name, rel_path = line.split("\t")
            frame = frame_map.get(orig_name)
            if frame:
                frame.path_mae = str((cwd / "prepared" / rel_path).resolve())
                frame.path_pdb = str((cwd / "prepared" / orig_name / "receptor.pdb").resolve())
                frame.state = "ready"
                prep_rec_count += 1

                system = session.query(System).filter_by(receptor_frame_id=frame.id).first()
                if system and system.pose:
                    prep_sdf = (cwd / "prepared" / orig_name / "native.sdf").resolve()
                    prep_mae = (cwd / "prepared" / orig_name / "native.mae").resolve()
                    if prep_sdf.exists() and prep_sdf.stat().st_size > 0:
                        system.pose.path_sdf = str(prep_sdf)
                        system.pose.path_mae = str(prep_mae) if prep_mae.exists() else None
                        system.pose.state = "ready"

    sub_receipt = cwd / "prep_ligand_receipt.tsv"
    if sub_receipt.exists():
        for line in sub_receipt.read_text().strip().splitlines():
            if not line.strip(): continue
            sub_name, rel_path = line.split("\t")
            sub_obj = session.query(Substance).filter_by(name=sub_name).first()
            if sub_obj:
                sub_obj.path_mae = str((cwd / "ligands" / rel_path).resolve())
                sdf_name = Path(rel_path).name.replace(".maegz", ".sdf").replace(".mae", ".sdf")
                sub_obj.path_sdf = str((cwd / "ligands" / Path(rel_path).parent / sdf_name).resolve())
                sub_obj.state = "ready"
                prep_sub_count += 1

    pose_receipt = cwd / "prep_pose_receipt.tsv"
    if pose_receipt.exists():
        for line in pose_receipt.read_text().strip().splitlines():
            if not line.strip(): continue
            pose_name, rel_mol2, rel_frcmod = line.split("\t")
            pose_obj = session.query(Pose).filter_by(name=pose_name).first()
            if pose_obj:
                pose_obj.path_mol2 = str((cwd / "prepared" / "prepared_poses" / rel_mol2).resolve())
                pose_obj.path_frcmod = str((cwd / "prepared" / "prepared_poses" / rel_frcmod).resolve())
                prep_pose_count += 1

    session.commit()
    session.close()

    console.print(f"\n[bold green]✓ Successfully prepared {prep_rec_count} receptors, {prep_sub_count} library substances, and parametrized {prep_pose_count} 3D poses![/bold green]")
