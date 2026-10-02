"""
File: cli/commands/dynamics.py
Description: CLI handler for 'sbdd md' using clean 3NF database schema.
Executes chunked MD and guarantees atomic registration of Trajectory records.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System, Trajectory
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_dynamics(args):
    cwd = Path.cwd()
    session = init_db()

    # 1. Resolve Systems needing simulation based on DB state
    reps_target = int(configs.get("md_simulation", {}).get("production", {}).get("replicates", 1))
    if getattr(args, "replicates", None):
        reps_target = int(args.replicates)

    query = session.query(System).filter_by(state="ready")
    if getattr(args, "system", None):
        query = query.filter_by(name=args.system)
    elif getattr(args, "target", None):
        target_obj = session.query(Target).filter_by(name=args.target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{args.target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter_by(target_id=target_obj.id)

    candidate_systems = query.all()

    # Skip systems that already have the required replicates completed
    systems = []
    for s in candidate_systems:
        existing_reps = session.query(Trajectory).filter_by(system_id=s.id).count()
        if existing_reps >= reps_target:
            continue
        systems.append(s)

    if not systems:
        if candidate_systems:
            console.print(f"[bold green]✓ All {len(candidate_systems)} ready System(s) already have {reps_target} replicates completed. Nothing to do![/bold green]")
        else:
            console.print("[yellow]No ready Systems found to simulate. Run 'sbdd prep' first.[/yellow]")
        session.close()
        return

    plan_list = []
    for sys_obj in systems:
        if sys_obj.state != "ready":
            console.print(f"[bold red]FATAL: System '{sys_obj.name}' is in state '{sys_obj.state}'. Must be 'ready'. Run 'sbdd prep'![/bold red]")
            sys.exit(1)

        frame = sys_obj.receptor_frame
        rec_path = Path(frame.path).resolve()

        lig_path = "NO_LIGAND"
        lig_name = ""
        if sys_obj.pose and Path(sys_obj.pose.path).exists():
            lig_path = str(Path(sys_obj.pose.path).resolve())
            lig_name = sys_obj.pose.ligand.name[:4].upper()

        plan_list.append(
            {
                "id": sys_obj.name,
                "rec_path": str(rec_path),
                "lig_path": lig_path,
                "lig_name": lig_name,
                "cofactors": frame.cofactors_list,
            }
        )

    # 2. Check MD protocol
    protocol_dir = cwd / "configs" / "md_protocol"
    if not protocol_dir.exists() or not (protocol_dir / "protocol.txt").exists():
        console.print("[bold cyan]MD protocol playlist not found. Auto-initializing 'sbdd md-init'...[/bold cyan]")
        from cli.main import cmd_md_init
        cmd_md_init(args)

    plan_file = cwd / "md_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(f"[bold cyan]Dispatched MD for {len(plan_list)} System(s): {[s['id'] for s in plan_list]}[/bold cyan]")

    # 3. Nextflow execution
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

    # 4. Atomic Registration of Trajectories into Database
    reps = int(configs.get("md_simulation", {}).get("production", {}).get("replicates", 1))
    if getattr(args, "replicates", None):
        reps = int(args.replicates)

    total_duration_ns = float(configs.get("md_simulation", {}).get("production", {}).get("total_time_ns", 0.0))
    if getattr(args, "total_ns", None):
        total_duration_ns = float(args.total_ns)

    registered_trajs = 0
    for sys_obj in systems:
        for rep_idx in range(1, reps + 1):
            traj_name = f"{sys_obj.name}_rep{rep_idx}"
            raw_dir = cwd / "md" / sys_obj.name / f"rep_{rep_idx}" / "raw"
            tpr_path = raw_dir / "sys_ref.tpr"

            # Check if simulation produced outputs
            chunks = list(raw_dir.glob("md_*.xtc"))
            if chunks and tpr_path.exists():
                ex_traj = session.query(Trajectory).filter_by(name=traj_name).first()
                if ex_traj:
                    ex_traj.system_id = sys_obj.id
                    ex_traj.replicate = rep_idx
                    ex_traj.length_ns = total_duration_ns
                    ex_traj.raw_dir = str(raw_dir.resolve())
                    ex_traj.path_tpr = str(tpr_path.resolve())
                    ex_traj.status = "raw"
                else:
                    new_traj = Trajectory(
                        system_id=sys_obj.id,
                        name=traj_name,
                        replicate=rep_idx,
                        length_ns=total_duration_ns,
                        raw_dir=str(raw_dir.resolve()),
                        path_tpr=str(tpr_path.resolve()),
                        status="raw",
                    )
                    session.add(new_traj)
                registered_trajs += 1

    session.commit()
    session.close()

    console.print(f"\n[bold green]✓ MD simulation completed! Registered {registered_trajs} Trajectory record(s) in DB.[/bold green]")
    console.print("[dim]Run 'sbdd traj-clean' to strip solvent and fix PBC for downstream analysis.[/dim]")
