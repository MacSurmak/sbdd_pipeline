"""
File: cli/commands/dynamics.py
Description: CLI handler for 'sbdd md' using clean 3NF database schema.
Executes chunked MD and guarantees atomic registration of Trajectory records.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Pose, System, Trajectory
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def run_dynamics(
    system: str | None = None,
    target: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """Executes chunked MD and guarantees atomic registration of Trajectory records."""
    cwd = Path.cwd()
    session = init_db()
    configs = load_project_configs(cwd)

    # 1. Resolve Systems needing simulation and calculate required top-up replicates
    # Strict Fail-Fast on missing parameters
    try:
        reps_target = int(configs["md_simulation"]["production"]["replicates"])
        total_duration_ns = float(
            configs["md_simulation"]["production"]["total_time_ns"]
        )
    except KeyError as e:
        console.print(
            f"[bold red]FATAL: Missing required MD configuration parameter in md.yaml: {e}[/bold red]"
        )
        sys.exit(1)

    query = session.query(System).filter_by(state="ready")
    if system:
        query = query.filter_by(name=system)
    elif target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter_by(target_id=target_obj.id)

    candidate_systems = query.all()

    # Determine systems with missing active replicates (active = raw or cleaned)
    plan_list = []
    systems_to_run = []

    for s in candidate_systems:
        existing_trajs = session.query(Trajectory).filter_by(system_id=s.id).all()

        # Purge any corrupted ghost records where raw_dir path doesn't match replicate number
        for t in existing_trajs:
            if f"rep_{t.replicate}" not in t.raw_dir:
                console.print(
                    f"[yellow]Purging corrupted trajectory record '{t.name}' (points to mismatched path {t.raw_dir})...[/yellow]"
                )
                session.delete(t)
                session.flush()

        # Identify truly valid existing replicates
        valid_reps = {
            t.replicate
            for t in existing_trajs
            if t.status in ["raw", "cleaned"] and f"rep_{t.replicate}" in t.raw_dir
        }

        # Calculate exact non-colliding replicate IDs needed to reach target
        missing_reps = [r for r in range(1, reps_target + 1) if r not in valid_reps]
        if not missing_reps:
            continue

        # Strict 3NF Fail-Fast: Dynamics requires fully parametrized topology
        if not s.path_gro or not s.path_top or not s.path_ndx:
            raise ValueError(
                f"FATAL: System '{s.name}' does not have prepared topology paths registered in 3NF database! "
                f"Entity must be parametrized cleanly via 'sbdd md-setup' before launching dynamics."
            )

        p_gro = Path(s.path_gro).resolve()
        p_top = Path(s.path_top).resolve()
        p_ndx = Path(s.path_ndx).resolve()

        for p in (p_gro, p_top, p_ndx):
            if not p.exists():
                raise FileNotFoundError(f"FATAL: Topology file missing on disk: {p}")

        # Strict Single Source of Truth: resolve all restraint files directly from DB
        env = s.environment_dict
        posre_files = env.get("posre_files", [])
        if not posre_files:
            raise ValueError(
                f"FATAL: No position restraint files found in DB for system '{s.name}'. "
                f"Was it setup with an older pipeline version?"
            )
            
        for f in posre_files:
            if not Path(f).exists():
                raise FileNotFoundError(f"FATAL: Registered position restraint missing on disk: {f}")

        systems_to_run.append((s, missing_reps))
        plan_list.append(
            {
                "id": s.name,
                "gro": str(p_gro),
                "top": str(p_top),
                "ndx": str(p_ndx),
                "posre_files": posre_files,
                "reps_to_run": missing_reps,
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

    total_missing_dispatched = sum(len(item["reps_to_run"]) for item in plan_list)
    max_reps_dispatch = max(len(item["reps_to_run"]) for item in plan_list)
    console.print(
        f"[bold cyan]Dispatching simulation to replenish replicates (Target: {reps_target}, Missing: {total_missing_dispatched})...[/bold cyan]"
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

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    ret_code = execute_nextflow(
        "dynamics", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow MD simulation failed.[/bold red]")
        sys.exit(ret_code)

    # 4. Atomic Registration of Trajectories directly linked to their published replicate paths
    registered_trajs = 0
    for sys_obj, missing_reps in systems_to_run:
        for rep_idx in missing_reps:
            raw_dir = cwd / "md_runs" / sys_obj.name / f"rep_{rep_idx}" / "raw"
            tpr_path = raw_dir / "sys_ref.tpr"
            first_chunk = raw_dir / "md_1.xtc"

            if not first_chunk.exists() or not tpr_path.exists():
                session.close()
                raise FileNotFoundError(
                    f"FATAL: Deterministic trajectory artifacts (md_1.xtc or sys_ref.tpr) missing for {sys_obj.name} replicate {rep_idx} in {raw_dir}. "
                    "Nextflow finished without generating expected deliverables."
                )

            traj_name = f"{sys_obj.name}_rep{rep_idx}"
            ex_t = session.query(Trajectory).filter_by(name=traj_name).first()
            if ex_t:
                ex_t.replicate = rep_idx
                ex_t.length_ns = total_duration_ns
                ex_t.raw_dir = str(raw_dir.resolve())
                ex_t.path_tpr = str(tpr_path.resolve())
                ex_t.status = "raw"
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

    console.print(
        f"\n[bold green]✓ MD simulation completed! Registered {registered_trajs} Trajectory record(s) in DB.[/bold green]"
    )
    console.print(
        "[dim]Run 'sbdd traj-clean' to strip solvent and fix PBC for downstream analysis.[/dim]"
    )
