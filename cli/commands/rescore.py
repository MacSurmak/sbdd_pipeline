"""
File: cli/commands/rescore.py
Description: CLI handler for 'sbdd rescore': standalone Prime MM-GBSA and Consensus rescoring.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from rich.table import Table
from cli.core.db import init_db, Target, ReceptorFrame, Pose
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_rescore(args):
    cwd = Path.cwd()
    session = init_db()

    dock_dir = cwd / "docking"
    if not dock_dir.exists():
        console.print(
            "[bold red]FATAL: No 'docking/' directory found. Run 'sbdd dock' first![/bold red]"
        )
        sys.exit(1)

    pv_files = sorted(dock_dir.glob("*/*_dock_pv.maegz"))
    if not pv_files:
        pv_files = sorted(dock_dir.glob("*_dock_pv.maegz"))

    if not pv_files:
        console.print(
            "[bold red]FATAL: No Glide Pose Viewer (*_dock_pv.maegz) files found to rescore.[/bold red]"
        )
        sys.exit(1)

    tasks = []
    for pv in pv_files:
        frame_name = pv.stem.replace("_dock_pv", "")
        tasks.append(
            {"name": frame_name, "pv_path": str(pv.resolve()), "frame": frame_name}
        )

    plan = {
        "use_whitelist": getattr(args, "use_whitelist", False),
        "whitelist_path": (
            str((cwd / "funnel_whitelist.csv").resolve())
            if (cwd / "funnel_whitelist.csv").exists()
            else ""
        ),
        "tasks": tasks,
    }

    plan_file = cwd / "rescore_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan, f, indent=2)

    console.print(
        f"[bold cyan]Dispatching standalone Prime MM-GBSA rescoring for {len(tasks)} receptor frame(s)...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="rescore")
    nf_args.extend(["--rescore_plan", str(plan_file.resolve())])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "rescore", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow MM-GBSA rescoring failed.[/bold red]")
        sys.exit(ret_code)

    console.print(
        "\n[bold green]✓ Prime MM-GBSA rescoring completed! Master ranking updated in reports/ranking.csv.[/bold green]\n"
    )
    session.close()
