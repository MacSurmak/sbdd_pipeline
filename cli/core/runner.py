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

def clean_shared_scratch(project_name: str, verbose: bool = True):
    """Safely cleans the project's scratch workspace on /shared."""
    scratch_dir = Path(f"/shared/work/{project_name}")
    if scratch_dir.exists():
        if verbose:
            console.print(f"[yellow]Cleaning temporary files in {scratch_dir}...[/yellow]")
        try:
            shutil.rmtree(scratch_dir)
            if verbose:
                console.print(f"[green]Successfully cleaned {scratch_dir}[/green]")
        except Exception as e:
            console.print(f"[red]Error cleaning {scratch_dir}: {e}[/red]")
    else:
        if verbose:
            console.print(f"[dim]No scratch directory found at {scratch_dir}[/dim]")

def execute_nextflow(entrypoint: str, extra_args: list, project_dir: Path,
                     profile: str = "slurm", resume: bool = True,
                     auto_clean: bool = True) -> int:
    """Invokes Nextflow orchestrator for a given subworkflow."""
    main_nf = FRAMEWORK_DIR / "main.nf"
    
    project_name = project_dir.name
    # Extract project_name from extra_args if present
    if "--project_name" in extra_args:
        idx = extra_args.index("--project_name")
        if idx + 1 < len(extra_args):
            project_name = extra_args[idx + 1]

    log_file = project_dir / "logs" / f"nextflow_{entrypoint}.log"
    (project_dir / "logs").mkdir(parents=True, exist_ok=True)

    cmd = [
        "nextflow",
        "-log", str(log_file.resolve()),
        "run", str(main_nf.resolve()),
        "-profile", profile,
        "--entry", entrypoint
    ]

    if resume:
        cmd.append("-resume")

    cmd.extend(extra_args)

    console.rule(f"[bold cyan]Nextflow Orchestrator: {entrypoint.upper()} ({profile.upper()})[/bold cyan]")
    console.print(f"[dim]Project Scratch: /shared/work/{project_name}[/dim]")
    console.print(f"[dim]Deliverables:    {project_dir.resolve()}[/dim]")
    console.print(f"[bold blue]Command:[/bold blue] {' '.join(cmd)}\n")

    try:
        proc = subprocess.Popen(cmd, cwd=project_dir)
        ret = proc.wait()

        if ret == 0:
            console.print(f"\n[bold green]✓ Subworkflow '{entrypoint}' finished successfully![/bold green]")
            if auto_clean:
                clean_shared_scratch(project_name, verbose=True)
            return 0
        else:
            console.print(f"\n[bold red]✗ Nextflow failed with exit code {ret}[/bold red]")
            console.print(f"[yellow]Scratch directory /shared/work/{project_name} preserved for '-resume'.[/yellow]")
            return ret

    except KeyboardInterrupt:
        console.print("\n[bold yellow]Execution canceled by user. Scratch preserved for resume.[/bold yellow]")
        return 130
    except Exception as e:
        console.print(f"\n[bold red]Failed to launch Nextflow: {e}[/bold red]")
        return 1

