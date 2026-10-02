"""
File: cli/commands/dynamics.py
Description: CLI handler for 'sbdd md' / 'sbdd dynamics'. Dispatches MD setup and chunked simulation on Slurm.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, System, Receptor, Pose, Ligand, Trajectory
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_dynamics(args):
    cwd = Path.cwd()
    session = init_db()

    # 1. Determine target System
    target_sys_name = getattr(args, "system", None)
    if target_sys_name:
        systems = session.query(System).filter_by(name=target_sys_name).all()
        if not systems:
            systems = session.query(System).filter_by(receptor_name=target_sys_name).all()
    else:
        # Default: all ready systems
        systems = session.query(System).filter_by(status="ready").all()

    if not systems:
        console.print(
            "[bold red]FATAL: No ready System found to simulate! Specify --system <name> from 'sbdd tree' or run 'sbdd prep'.[/bold red]"
        )
        sys.exit(1)

    plan_list = []
    for sys_obj in systems:
        rec = session.query(Receptor).filter_by(name=sys_obj.receptor_name).first()
        if not rec or rec.status != "ready":
            console.print(
                f"[bold red]FATAL: Receptor '{sys_obj.receptor_name}' for system '{sys_obj.name}' is not ready! Run 'sbdd prep' first.[/bold red]"
            )
            sys.exit(1)

        rec_path = Path(rec.path).resolve()

        # Ligand handling
        lig_path = "NO_LIGAND"
        lig_name = ""
        if sys_obj.pose_name:
            pose = session.query(Pose).filter_by(name=sys_obj.pose_name).first()
            if pose and Path(pose.path).exists():
                lig_path = str(Path(pose.path).resolve())
                lig_name = pose.ligand.meta_dict.get("resname", pose.ligand_name)[:4].upper() if pose.ligand else "LIG"

        cofactors = rec.meta_dict.get("cofactors", [])

        plan_list.append(
            {
                "id": sys_obj.name,
                "rec_path": str(rec_path),
                "lig_path": lig_path,
                "lig_name": lig_name,
                "cofactors": cofactors,
            }
        )

    # 2. Check MD protocol directory
    protocol_dir = cwd / "configs" / "md_protocol"
    if not protocol_dir.exists() or not (protocol_dir / "protocol.txt").exists():
        console.print(
            "[bold cyan]MD protocol playlist not found. Auto-initializing 'sbdd md-init'...[/bold cyan]"
        )
        from cli.main import cmd_md_init
        cmd_md_init(args)

    plan_file = cwd / "md_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Generated md_plan.json for {len(plan_list)} System(s): {[p['id'] for p in plan_list]}[/bold cyan]"
    )

    # 3. Build Nextflow args
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="dynamics")
    nf_args.extend(["--md_plan", str(plan_file.resolve())])

    if getattr(args, "replicates", None):
        if "--md_replicates" in nf_args:
            idx = nf_args.index("--md_replicates")
            nf_args[idx + 1] = str(args.replicates)
        else:
            nf_args.extend(["--md_replicates", str(args.replicates)])

    if getattr(args, "total_ns", None):
        if "--md_total_ns" in nf_args:
            idx = nf_args.index("--md_total_ns")
            nf_args[idx + 1] = str(args.total_ns)
        else:
            nf_args.extend(["--md_total_ns", str(args.total_ns)])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "dynamics", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow MD simulation failed.[/bold red]")
        sys.exit(ret_code)

    # 4. Automatically register completed trajectories into the database
    reps = int(configs.get("md_simulation", {}).get("production", {}).get("replicates", 1))
    if getattr(args, "replicates", None):
        reps = int(args.replicates)

    total_duration_ns = float(configs.get("md_simulation", {}).get("production", {}).get("total_time_ns", 0.0))
    if getattr(args, "total_ns", None):
        total_duration_ns = float(args.total_ns)

    registered_trajs = 0
    for plan_item in plan_list:
        sys_name = plan_item["id"]
        rec_id = sys_name.replace("_ready_holo", "").replace("_holo", "").replace("_ready_apo", "").replace("_apo", "")

        for rep_idx in range(1, reps + 1):
            traj_name = f"{rec_id}_rep{rep_idx}"
            raw_dir = cwd / "md" / sys_name / f"rep_{rep_idx}" / "raw"
            if not raw_dir.exists():
                raw_dir = cwd / "md" / rec_id / f"rep_{rep_idx}" / "raw"

            tpr_path = raw_dir / "sys_ref.tpr"

            chunks = sorted(
                raw_dir.glob("md_*.xtc"),
                key=lambda p: int("".join(filter(str.isdigit, p.stem)) or 0)
            )

            if chunks and tpr_path.exists():
                primary_xtc = chunks[0]
                ex_traj = session.query(Trajectory).filter_by(name=traj_name).first()
                if ex_traj:
                    ex_traj.system_name = sys_name
                    ex_traj.replicate = rep_idx
                    ex_traj.path_xtc = str(primary_xtc.resolve())
                    ex_traj.path_tpr = str(tpr_path.resolve())
                    ex_traj.length_ns = total_duration_ns
                    ex_traj.status = "raw"
                else:
                    new_traj = Trajectory(
                        name=traj_name,
                        system_name=sys_name,
                        replicate=rep_idx,
                        path_xtc=str(primary_xtc.resolve()),
                        path_tpr=str(tpr_path.resolve()),
                        length_ns=total_duration_ns,
                        status="raw"
                    )
                    session.add(new_traj)
                registered_trajs += 1

    session.commit()
    console.print(
        f"\n[bold green]✓ MD simulation finished successfully! Registered {registered_trajs} Trajectory record(s) in DB.[/bold green]"
    )
    console.print("[dim]Use 'sbdd traj-clean' to strip solvent and fix PBC for downstream analysis.[/dim]")
    session.close()
