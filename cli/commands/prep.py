"""
File: cli/commands/prep.py
Description: CLI handler for 'sbdd prep' using clean 3NF database lineage.
Prepares raw ReceptorFrames and in-situ Poses into ready states without mutating raw records.
"""

import json
from pathlib import Path
import sys
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_prep(args):
    cwd = Path.cwd()
    session = init_db()

    # 1. Fetch raw frames that need preparation
    query = session.query(ReceptorFrame).filter(ReceptorFrame.state == "raw")
    if getattr(args, "target", None):
        target_obj = session.query(Target).filter_by(name=args.target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{args.target}' not found in database.[/bold red]")
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif getattr(args, "frames", None):
        frame_names = [f.strip() for f in args.frames.split(",") if f.strip()]
        query = query.filter(ReceptorFrame.name.in_(frame_names))

    raw_frames = query.all()
    if not raw_frames:
        console.print("[green]No raw frames needing preparation found. Everything is up to date![/green]")
        session.close()
        return

    # 2. Build Plan for Nextflow
    plan_data = []
    frame_map = {}

    for frame in raw_frames:
        frame_map[frame.name] = frame

        # Find in-situ pose in holo system for this frame if any
        system = session.query(System).filter_by(receptor_frame_id=frame.id).first()
        native_lig_code = ""
        if system and system.pose:
            native_lig_code = system.pose.ligand.name

        prep_mode = "light" if frame.source_type == "md_centroid" else "full"

        plan_data.append(
            {
                "name": frame.name,
                "path": str(Path(frame.path).resolve()),
                "prep_mode": prep_mode,
                "cofactors": frame.cofactors_list,
                "native_ligand": native_lig_code,
            }
        )

    plan_file = cwd / "prep_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_data, f, indent=2)

    console.print(f"[bold cyan]Generated prep_plan.json for {len(raw_frames)} frame(s)...[/bold cyan]")

    # 3. Build Nextflow args and execute
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="prepare")
    nf_args.extend(["--plan", str(plan_file.resolve())])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "prepare", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow preparation failed. Database not updated.[/bold red]")
        sys.exit(ret_code)

    # 4. Read Receipt and instantiate ready records with data lineage
    receipt_file = cwd / "prep_receipt.tsv"
    if not receipt_file.exists():
        console.print("[bold red]FATAL: prep_receipt.tsv not found! Nextflow execution might have failed.[/bold red]")
        sys.exit(1)

    prepared_count = 0
    lines = receipt_file.read_text().strip().splitlines()

    for line in lines:
        if not line.strip():
            continue
        orig_frame_name, raw_out_path = line.split("\t")
        raw_frame = frame_map.get(orig_frame_name)
        if not raw_frame:
            continue

        prep_frame_path = Path(raw_out_path)
        if not prep_frame_path.exists():
            # Check published location
            cand = cwd / "prepared" / f"{orig_frame_name}_prepared.mae"
            if cand.exists():
                prep_frame_path = cand

        # A. Register new ready ReceptorFrame derived from raw_frame
        prep_frame_name = f"{orig_frame_name}_prep"
        ready_frame = session.query(ReceptorFrame).filter_by(name=prep_frame_name).first()
        if not ready_frame:
            ready_frame = ReceptorFrame(
                target_id=raw_frame.target_id,
                name=prep_frame_name,
                source_type=raw_frame.source_type,
                state="ready",
                parent_frame_id=raw_frame.id,
                path=str(prep_frame_path.resolve()),
            )
            ready_frame.cofactors_list = raw_frame.cofactors_list
            session.add(ready_frame)
            session.flush()
        else:
            ready_frame.path = str(prep_frame_path.resolve())
            ready_frame.state = "ready"
            ready_frame.cofactors_list = raw_frame.cofactors_list

        # B. Check for prepared in-situ native ligand
        prep_pose_obj = None
        raw_system = session.query(System).filter_by(receptor_frame_id=raw_frame.id).first()

        if raw_system and raw_system.pose:
            raw_pose = raw_system.pose
            prep_sdf_cand = cwd / "prepared" / f"{orig_frame_name}_ref_native_prepared.sdf"
            if not prep_sdf_cand.exists():
                prep_sdf_cand = cwd / "prepared" / f"{orig_frame_name}_ref_native.sdf"

            if prep_sdf_cand.exists() and prep_sdf_cand.stat().st_size > 0:
                prep_pose_name = f"{raw_pose.name}_prep"
                ready_pose = session.query(Pose).filter_by(name=prep_pose_name).first()
                if not ready_pose:
                    ready_pose = Pose(
                        ligand_id=raw_pose.ligand_id,
                        receptor_frame_id=ready_frame.id,
                        name=prep_pose_name,
                        origin="crystal",
                        state="ready",
                        parent_pose_id=raw_pose.id,
                        path=str(prep_sdf_cand.resolve()),
                    )
                    session.add(ready_pose)
                    session.flush()
                else:
                    ready_pose.path = str(prep_sdf_cand.resolve())
                    ready_pose.state = "ready"
                prep_pose_obj = ready_pose

        # C. Create Ready System (Apo or Holo)
        ready_sys_name = f"{orig_frame_name}_holo" if prep_pose_obj else f"{orig_frame_name}_apo"
        ready_sys = session.query(System).filter_by(name=ready_sys_name).first()
        if not ready_sys:
            ready_sys = System(
                target_id=raw_frame.target_id,
                receptor_frame_id=ready_frame.id,
                pose_id=prep_pose_obj.id if prep_pose_obj else None,
                name=ready_sys_name,
                state="ready",
                origin=raw_frame.source_type,
            )
            session.add(ready_sys)
        else:
            ready_sys.receptor_frame_id = ready_frame.id
            ready_sys.pose_id = prep_pose_obj.id if prep_pose_obj else None
            ready_sys.state = "ready"

        prepared_count += 1

    session.commit()
    session.close()

    console.print(f"\n[bold green]✓ Successfully prepared {prepared_count} frame(s) into ready systems![/bold green]")
    console.print("[dim]Original raw frames and crystal poses preserved. Use 'sbdd tree' to view lineage.[/dim]")
