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
    configs = load_project_configs(cwd)

    # 1. Resolve Systems needing simulation and calculate required top-up replicates
    # Strict Fail-Fast on missing parameters
    try:
        reps_target = int(configs["md_simulation"]["production"]["replicates"])
        total_duration_ns = float(configs["md_simulation"]["production"]["total_time_ns"])
    except KeyError as e:
        console.print(f"[bold red]FATAL: Missing required MD configuration parameter in md.yaml: {e}[/bold red]")
        sys.exit(1)

    query = session.query(System).filter_by(state="ready")
    if getattr(args, "system", None):
        query = query.filter_by(name=args.system)
    elif getattr(args, "target", None):
        target_obj = session.query(Target).filter_by(name=args.target).first()
        if not target_obj:
            console.print(
                f"[bold red]FATAL: Target '{args.target}' not found.[/bold red]"
            )
            sys.exit(1)
        query = query.filter_by(target_id=target_obj.id)

    candidate_systems = query.all()

    # Determine systems with missing active replicates (active = raw or cleaned)
    plan_list = []
    systems_to_run = []

    for s in candidate_systems:
        active_reps = (
            session.query(Trajectory)
            .filter(
                Trajectory.system_id == s.id, Trajectory.status.in_(["raw", "cleaned"])
            )
            .count()
        )

        needed_reps = max(0, reps_target - active_reps)
        if needed_reps <= 0:
            continue

        frame = s.receptor_frame
        rec_path = Path(frame.path).resolve()

        lig_path = "NO_LIGAND"
        lig_name = ""
        if s.pose and Path(s.pose.path).exists():
            lig_path = str(Path(s.pose.path).resolve())
            lig_name = s.pose.ligand.name[:4].upper()

        systems_to_run.append((s, needed_reps))
        plan_list.append(
            {
                "id": s.name,
                "rec_path": str(rec_path),
                "lig_path": lig_path,
                "lig_name": lig_name,
                "cofactors": frame.cofactors_list,
                "reps_to_run": needed_reps,
            }
        )

    if not plan_list:
        if candidate_systems:
            console.print(
                f"[bold green]✓ All {len(candidate_systems)} ready System(s) already have {reps_target} active replicates. Nothing to do![/bold green]"
            )
        else:
            console.print(
                "[yellow]No ready Systems found to simulate. Run 'sbdd prep' first.[/yellow]"
            )
        session.close()
        return

    # Maximum reps needed across queued systems
    max_reps_dispatch = max(item["reps_to_run"] for item in plan_list)
    console.print(
        f"[bold cyan]Dispatching simulation to replenish replicates (Target: {reps_target}, Max to run: {max_reps_dispatch})...[/bold cyan]"
    )

    # 2. Check MD protocol
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
        f"[bold cyan]Dispatched MD for {len(plan_list)} System(s): {[s['id'] for s in plan_list]}[/bold cyan]"
    )

    # 3. Nextflow execution
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="dynamics")
    nf_args.extend(["--md_plan", str(plan_file.resolve())])

    # Enforce dynamic maximum replicates needed for this specific batch execution
    if "--md_replicates" in nf_args:
        idx = nf_args.index("--md_replicates")
        nf_args[idx + 1] = str(max_reps_dispatch)
    else:
        nf_args.extend(["--md_replicates", str(max_reps_dispatch)])

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

    # 4. Atomic Registration of Trajectories with Incremental Replicate Indexing
    registered_trajs = 0
    for sys_obj, needed_reps in systems_to_run:
        # Find highest existing replicate index
        existing_trajs = session.query(Trajectory).filter_by(system_id=sys_obj.id).all()
        max_idx = max([t.replicate for t in existing_trajs], default=0)

        for offset in range(1, needed_reps + 1):
            run_rep_idx = offset  # Index in current Nextflow run output
            new_global_rep_idx = max_idx + offset  # Non-colliding global index

            raw_dir = cwd / "md" / sys_obj.name / f"rep_{run_rep_idx}" / "raw"
            tpr_path = raw_dir / "sys_ref.tpr"
            chunks = list(raw_dir.glob("md_*.xtc"))

            if chunks and tpr_path.exists():
                traj_name = f"{sys_obj.name}_rep{new_global_rep_idx}"
                new_traj = Trajectory(
                    system_id=sys_obj.id,
                    name=traj_name,
                    replicate=new_global_rep_idx,
                    length_ns=total_duration_ns,
                    raw_dir=str(raw_dir.resolve()),
                    path_tpr=str(tpr_path.resolve()),
                    status="raw",
                )
                session.add(new_traj)
                registered_trajs += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ MD simulation completed! Registered {registered_trajs} Trajectory record(s) in DB.[/bold green]"
    )
    console.print(
        "[dim]Run 'sbdd traj-clean' to strip solvent and fix PBC for downstream analysis.[/dim]"
    )
