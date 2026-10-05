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
            console.print(
                f"[bold red]FATAL: Target '{args.target}' not found in database.[/bold red]"
            )
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif getattr(args, "frames", None):
        frame_names = [f.strip() for f in args.frames.split(",") if f.strip()]
        query = query.filter(ReceptorFrame.name.in_(frame_names))

    raw_frames = query.all()
    if not raw_frames:
        console.print(
            "[green]No raw frames needing preparation found. Everything is up to date![/green]"
        )
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

    console.print(
        f"[bold cyan]Generated prep_plan.json for {len(raw_frames)} frame(s)...[/bold cyan]"
    )

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
        console.print(
            "[bold red]Nextflow preparation failed. Database not updated.[/bold red]"
        )
        sys.exit(ret_code)

    # 4. Read Receipt and instantiate ready records with data lineage
    receipt_file = cwd / "prep_receipt.tsv"
    if not receipt_file.exists():
        console.print(
            "[bold red]FATAL: prep_receipt.tsv not found! Nextflow execution might have failed.[/bold red]"
        )
        sys.exit(1)

    prepared_count = 0
    lines = receipt_file.read_text().strip().splitlines()

    for line in lines:
        if not line.strip():
            continue
        orig_frame_name, mae_filename = line.split("\t")
        frame = frame_map.get(orig_frame_name)
        if not frame:
            continue

        # Strictly construct and enforce permanent project deliverable location
        prep_mae_path = (cwd / "prepared" / Path(mae_filename).name).resolve()
        if not prep_mae_path.exists():
            raise FileNotFoundError(
                f"FATAL: Prepared receptor MAE not found at canonical location: {prep_mae_path}"
            )

        frame.path = str(prep_mae_path)
        frame.state = "ready"

        # Update existing Pose in-place if native ligand was prepared
        system = session.query(System).filter_by(receptor_frame_id=frame.id).first()
        if system and system.pose:
            prep_sdf_path = (cwd / "prepared" / f"{orig_frame_name}_ref_native_prepared.sdf").resolve()
            if not prep_sdf_path.exists():
                prep_sdf_path = (cwd / "prepared" / f"{orig_frame_name}_ref_native.sdf").resolve()

            if prep_sdf_path.exists() and prep_sdf_path.stat().st_size > 0:
                system.pose.path = str(prep_sdf_path)
                system.pose.state = "ready"

        if system:
            system.state = "ready"

        prepared_count += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Successfully prepared {prepared_count} frame(s) into ready systems![/bold green]"
    )
    console.print(
        "[dim]Original raw frames and crystal poses preserved. Use 'sbdd tree' to view lineage.[/dim]"
    )
