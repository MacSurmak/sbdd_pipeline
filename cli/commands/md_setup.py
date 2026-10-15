"""
File: cli/commands/md_setup.py
Description: CLI handler for 'sbdd md-setup'.
Gathers 'raw' Systems, generates topologies, and atomically registers them as 'ready'.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Target, System
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def run_md_setup(
    system: str | None = None,
    target: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    cwd = Path.cwd()
    session = init_db()
    configs = load_project_configs(cwd)

    query = session.query(System).filter_by(state="raw")
    if system:
        query = query.filter_by(name=system)
    elif target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter_by(target_id=target_obj.id)

    raw_systems = query.all()
    if not raw_systems:
        console.print(
            "[green]No raw Systems found to setup. Topologies are up to date![/green]"
        )
        session.close()
        return

    plan_list = []
    for s in raw_systems:
        frame = s.receptor_frame
        if not frame.path_mae:
            raise ValueError(
                f"FATAL: ReceptorFrame '{frame.name}' does not have path_mae registered. "
                "Run 'sbdd prep' first."
            )
        rec_path = Path(frame.path_mae).resolve()
        if not rec_path.exists():
            raise FileNotFoundError(f"Prepared receptor missing: {rec_path}")

        lig_path = None
        lig_name = ""
        if s.pose:
            pose_sdf = s.pose.path_sdf or s.pose.path
            if not pose_sdf:
                raise ValueError(
                    f"FATAL: Native Pose '{s.pose.name}' missing path_sdf."
                )
            p_sdf_path = Path(pose_sdf).resolve()
            if not p_sdf_path.exists():
                raise FileNotFoundError(f"Prepared ligand missing: {p_sdf_path}")
            lig_path = str(p_sdf_path)
            lig_name = s.pose.ligand.name[:4].upper()

        plan_list.append(
            {
                "id": s.name,
                "rec_path": str(rec_path),
                "lig_path": lig_path,
                "lig_name": lig_name,
                "cofactors": frame.cofactors_list,
            }
        )

    plan_file = cwd / "md_setup_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Dispatching topology setup for {len(plan_list)} System(s)...[/bold cyan]"
    )

    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-setup")
    nf_args.extend(["--md_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
    ret_code = execute_nextflow(
        "md-setup", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow MD setup failed.[/bold red]")
        sys.exit(ret_code)

    receipt_file = cwd / "md_setup_receipt.tsv"
    if not receipt_file.exists():
        console.print("[bold red]FATAL: md_setup_receipt.tsv not found.[/bold red]")
        sys.exit(1)

    sys_map = {s.name: s for s in raw_systems}
    registered = 0
    lines = receipt_file.read_text().strip().splitlines()
    for line in lines:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 4:
            continue
        sys_id, gro, top, ndx = parts
        if sys_id in sys_map:
            s = sys_map[sys_id]
            s.path_gro = gro
            s.path_top = top
            s.path_ndx = ndx
            s.state = "ready"
            registered += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ MD Setup completed! {registered} System(s) transitioned to 'ready'.[/bold green]"
    )
    console.print(
        "[dim]Topologies compiled. Run 'sbdd md' to launch production dynamics.[/dim]"
    )
