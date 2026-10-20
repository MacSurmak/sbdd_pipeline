"""
File: cli/core/runner.py
Description: Nextflow orchestrator execution wrapper with Slurm monitoring and /shared scratch lifecycle.
"""

import subprocess
import shutil
import sys
from pathlib import Path
from rich.console import Console

console = Console()
FRAMEWORK_DIR = Path(__file__).resolve().parent.parent.parent


def clean_shared_scratch(workdir: str, verbose: bool = True):
    """Safely cleans the project's scratch workspace."""
    scratch_dir = Path(workdir)
    if scratch_dir.exists():
        if verbose:
            console.print(
                f"[yellow]Cleaning temporary files in {scratch_dir}...[/yellow]"
            )
        try:
            shutil.rmtree(scratch_dir)
            if verbose:
                console.print(f"[green]Successfully cleaned {scratch_dir}[/green]")
        except Exception as e:
            console.print(f"[red]Error cleaning {scratch_dir}: {e}[/red]")
    else:
        if verbose:
            console.print(f"[dim]No scratch directory found at {scratch_dir}[/dim]")


def execute_nextflow(
    entrypoint: str,
    extra_args: list,
    project_dir: Path,
    profile: str = "slurm",
    resume: bool = True,
    auto_clean: bool = False,
    configs: dict = None,
) -> int:
    """Invokes Nextflow orchestrator directly from the pristine /shared/sbdd_pipeline source tree."""
    import os

    target_main = FRAMEWORK_DIR / "main.nf"

    project_name = project_dir.name
    if "--project_name" in extra_args:
        idx = extra_args.index("--project_name")
        if idx + 1 < len(extra_args):
            project_name = extra_args[idx + 1]

    # Resolve scratch base safely via config or fallback
    default_scratch = f"/shared/work/{project_name}"
    scratch_base_str = default_scratch
    if configs:
        scratch_base_str = configs.get("project", {}).get("workdir", default_scratch)

    scratch_base = Path(scratch_base_str)
    scratch_nxf_home = scratch_base / ".nextflow"
    scratch_work = scratch_base / "work"
    scratch_nxf_home.mkdir(parents=True, exist_ok=True)
    scratch_work.mkdir(parents=True, exist_ok=True)

    reports_dir = project_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    log_file = reports_dir / f"nextflow_{entrypoint}.log"

    env = os.environ.copy()
    env["NXF_HOME"] = str(scratch_nxf_home.resolve())
    env["NXF_WORK"] = str(scratch_work.resolve())

    cmd = [
        "nextflow",
        "-log",
        str(log_file.resolve()),
        "run",
        str(target_main.resolve()),
        "-work-dir",
        str(scratch_work.resolve()),
        "-profile",
        profile,
        "--entry",
        entrypoint,
    ]

    if resume:
        cmd.append("-resume")

    cmd.extend(extra_args)

    console.rule(
        f"[bold cyan]Nextflow Orchestrator: {entrypoint.upper()} ({profile.upper()})[/bold cyan]"
    )
    console.print(f"[dim]Project Scratch: {scratch_base}[/dim]")
    console.print(f"[dim]Deliverables:    {project_dir.resolve()}[/dim]")
    console.print(f"[bold blue]Command:[/bold blue] {' '.join(cmd)}\n")

    # Run Nextflow INSIDE the scratch directory.
    # This prevents .nextflow/ caches and .nextflow.log from polluting the user's clean project folder.
    try:
        proc = subprocess.Popen(cmd, cwd=scratch_base, env=env)
        ret = proc.wait()

        if ret == 0:
            console.print(
                f"\n[bold green]✓ Subworkflow '{entrypoint}' finished successfully![/bold green]"
            )
            if auto_clean:
                clean_shared_scratch(str(scratch_base), verbose=True)
            return 0
        else:
            console.print(
                f"\n[bold red]✗ Nextflow failed with exit code {ret}[/bold red]"
            )
            console.print(
                "[dim]Note: Tasks execute in node-local scratch ($SLURM_TMPDIR) and are purged upon exit. Intermediate calculation artifacts do not reside in /shared/work/.[/dim]"
            )
            console.print(
                f"[yellow]Scratch directory {scratch_base} preserved for '-resume'.[/yellow]"
            )
            return ret

    except KeyboardInterrupt:
        console.print(
            "\n[bold yellow]Execution canceled by user. Scratch preserved for resume.[/bold yellow]"
        )
        return 130
    except Exception as e:
        console.print(f"\n[bold red]Failed to launch Nextflow: {e}[/bold red]")
        return 1
