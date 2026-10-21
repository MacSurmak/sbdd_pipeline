"""
File: cli/commands/project.py
Description: Project initialization, monitoring, and cleaning utilities.
"""

import sys
import shutil
import subprocess
from pathlib import Path
from rich.console import Console
from cli.core.config_loader import load_project_configs
from cli.core.runner import clean_shared_scratch, FRAMEWORK_DIR

console = Console()


def run_init(
    directory: str = ".", name: str | None = None, force: bool = False
) -> None:
    import random
    import string

    target_dir = Path(directory).resolve()
    target_dir.mkdir(parents=True, exist_ok=True)

    config_dir = target_dir / "configs"
    config_dir.mkdir(parents=True, exist_ok=True)

    project_name = name
    if not project_name:
        suffix = "".join(random.choices(string.ascii_lowercase + string.digits, k=6))
        project_name = f"SBDD_{suffix}"

    src_configs = FRAMEWORK_DIR / "assets" / "configs"
    if src_configs.exists():
        for f in src_configs.glob("*.yaml"):
            dest = config_dir / f.name
            if not dest.exists() or force:
                if f.name == "project.yaml":
                    content = f.read_text(encoding="utf-8")
                    content = content.replace(
                        'name: "SBDD_Project"', f'name: "{project_name}"'
                    )
                    dest.write_text(content, encoding="utf-8")
                else:
                    shutil.copy(f, dest)

    from cli.core.db import init_db

    init_db(target_dir / "sbdd_state.db", create_if_missing=True)

    console.print(
        f"[bold green]✓ Project '{project_name}' initialized at: {target_dir}[/bold green]"
    )
    console.print(f"[cyan]Configurations created in:[/cyan] {config_dir}")
    console.print("\n[bold]Next steps:[/bold]")
    console.print(
        "  1. Copy your structures (.cif or .pdb) into the project directory."
    )
    console.print("  2. Edit configs/manifest.yaml to describe your targets.")
    console.print("  3. Review configurations (e.g., md.yaml, gnina.yaml).")
    console.print(
        "  4. Run [bold cyan]sbdd prep[/bold cyan] to standardize the structures!\n"
    )


def run_status() -> None:
    console.rule("[bold cyan]Cluster & Slurm Status[/bold cyan]")
    try:
        res = subprocess.run(
            ["squeue", "-o", "%.10i %.10P %.20j %.10u %.2t %.10M %.6D %R"],
            capture_output=True,
            text=True,
        )
        console.print(res.stdout)
    except Exception as e:
        console.print(f"[red]Could not query Slurm squeue: {e}[/red]")

    log_f = Path.cwd() / ".nextflow.log"
    if log_f.exists():
        console.rule("[bold cyan]Recent Nextflow Activity[/bold cyan]")
        tail = subprocess.run(
            ["tail", "-n", "15", str(log_f)], capture_output=True, text=True
        )
        console.print(tail.stdout)


def run_md_init() -> None:
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    md_dir = cwd / "configs" / "md_protocol"
    md_dir.mkdir(parents=True, exist_ok=True)

    is_mem = "true" if configs.get("membrane", {}).get("enabled", False) else "false"
    tpl_dir = FRAMEWORK_DIR / "assets" / "mdp"
    script = FRAMEWORK_DIR / "bin" / "generate_mdp.py"

    temp = configs.get("md_simulation", {}).get("production", {}).get(
        "temperature"
    ) or configs.get("mmpbsa", {}).get("temperature")
    if not temp:
        console.print(
            "[bold red]FATAL: Missing simulation temperature in configs (md.yaml or mmpbsa.yaml)![/bold red]"
        )
        sys.exit(1)

    prod = configs.get("md_simulation", {}).get("production", {})
    chunk_ns = prod.get("chunk_size_ns")
    if not chunk_ns:
        console.print(
            "[bold red]FATAL: Missing 'md_simulation.production.chunk_size_ns' in configs/md.yaml![/bold red]"
        )
        sys.exit(1)

    cmd = [
        sys.executable,
        str(script),
        "--templates_dir",
        str(tpl_dir),
        "--outdir",
        str(md_dir),
        "--is_membrane",
        is_mem,
        "--temperature",
        str(temp),
        "--chunk_ns",
        str(chunk_ns),
        "--dt_ps",
        str(prod.get("dt_ps", "0.002")),
        "--compressibility",
        str(prod.get("compressibility", "4.5e-5")),
        "--tau_p",
        str(prod.get("tau_p", "5.0")),
        "--prod_restraints_fc",
        str(prod.get("restraints_fc", "")),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        console.print(
            f"[bold red]Failed to generate MD protocol: {res.stderr}[/bold red]"
        )
        sys.exit(1)

    console.print(
        f"[bold green]✓ MD Protocol initialized from assets/mdp into {md_dir}[/bold green]"
    )
    console.print(
        "[dim]You can now edit the generated .mdp files or reorder protocol.txt before running MD.[/dim]"
    )


def run_clean() -> None:
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    proj_name = configs.get("project", {}).get("name", cwd.name)
    workdir = configs.get("project", {}).get("workdir", f"/shared/work/{proj_name}")
    clean_shared_scratch(workdir, verbose=True)
