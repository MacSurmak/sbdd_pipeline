"""
File: cli/main.py
Description: Main CLI entrypoint for `sbdd` orchestrator tool.
"""
import argparse
import sys
import shutil
import subprocess
from pathlib import Path
from rich.console import Console
from rich.table import Table

from cli.core.config_loader import load_project_configs, build_nextflow_args
from cli.core.runner import execute_nextflow, clean_shared_scratch, FRAMEWORK_DIR

console = Console()

def cmd_init(args):
    """Initializes a new project workspace with atomic YAML configs."""
    target_dir = Path(args.directory).resolve()
    target_dir.mkdir(parents=True, exist_ok=True)
    
    config_dir = target_dir / "configs"
    config_dir.mkdir(parents=True, exist_ok=True)

    src_configs = FRAMEWORK_DIR / "assets" / "configs"
    if src_configs.exists():
        for f in src_configs.glob("*.yaml"):
            dest = config_dir / f.name
            if not dest.exists() or args.force:
                shutil.copy(f, dest)

    console.print(f"[bold green]✓ Project initialized at: {target_dir}[/bold green]")
    console.print(f"[cyan]Configurations created in:[/cyan] {config_dir}")
    console.print("\n[bold]Next steps:[/bold]")
    console.print("  1. Copy your receptor (.cif or .pdb) into the project directory.")
    console.print("  2. Copy your ligand library (.sdf) into the project directory.")
    console.print("  3. Review configs in configs/ (e.g. receptor.yaml, gnina.yaml).")
    console.print("  4. Run [bold cyan]sbdd prep[/bold cyan] or [bold cyan]sbdd all[/bold cyan]!\n")

def cmd_status(args):
    """Inspects active Slurm jobs and Nextflow status."""
    console.rule("[bold cyan]Cluster & Slurm Status[/bold cyan]")
    try:
        res = subprocess.run(["squeue", "-o", "%.10i %.10P %.20j %.10u %.2t %.10M %.6D %R"], capture_output=True, text=True)
        console.print(res.stdout)
    except Exception as e:
        console.print(f"[red]Could not query Slurm squeue: {e}[/red]")

    # Check Nextflow log in cwd
    log_f = Path.cwd() / ".nextflow.log"
    if log_f.exists():
        console.rule("[bold cyan]Recent Nextflow Activity[/bold cyan]")
        tail = subprocess.run(["tail", "-n", "15", str(log_f)], capture_output=True, text=True)
        console.print(tail.stdout)

def cmd_clean(args):
    """Cleans /shared/work/<project_name> scratch directory."""
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    proj_name = configs.get("project", {}).get("name", cwd.name)
    clean_shared_scratch(proj_name, verbose=True)

def run_workflow(entrypoint: str, args):
    """Generic dispatcher from CLI flags and configs to Nextflow execution."""
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd)

    # CLI flag overrides
    if getattr(args, "receptor", None):
        nf_args.extend(["--receptor", str(Path(args.receptor).resolve())])
    if getattr(args, "library", None):
        nf_args.extend(["--library", str(Path(args.library).resolve())])
    if getattr(args, "engine", None):
        nf_args.extend(["--dock_engine", args.engine])
    if getattr(args, "hits", None):
        nf_args.extend(["--md_hits", args.hits])
    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)
    auto_clean = not getattr(args, "keep_scratch", False)

    return execute_nextflow(entrypoint, nf_args, cwd, profile=profile, resume=resume, auto_clean=auto_clean)

def main():
    parser = argparse.ArgumentParser(
        prog="sbdd",
        description="SBDD Pipeline Nextflow HPC Orchestrator CLI"
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # 1. init
    p_init = subparsers.add_parser("init", help="Initialize new project workspace")
    p_init.add_argument("directory", nargs="?", default=".", help="Target project directory (default: current)")
    p_init.add_argument("--force", action="store_true", help="Overwrite existing templates")

    # Common options for workflows
    def add_common_wf_args(p):
        p.add_argument("--local", action="store_true", help="Run locally without Slurm")
        p.add_argument("--no-resume", action="store_true", help="Start from scratch without -resume")
        p.add_argument("--keep-scratch", action="store_true", help="Do not delete /shared/work upon completion")
        p.add_argument("--urgent", action="store_true", help="Route tasks to urgent Slurm partition")

    # 2. prep
    p_prep = subparsers.add_parser("prep", help="Run receptor & ligand preparation")
    p_prep.add_argument("--receptor", help="Path to input structure (.cif/.pdb)")
    p_prep.add_argument("--library", help="Path to ligand library (.sdf)")
    add_common_wf_args(p_prep)

    # 3. dock
    p_dock = subparsers.add_parser("dock", help="Run distributed docking across GPU nodes")
    p_dock.add_argument("--engine", choices=["gnina", "glide"], help="Docking engine choice")
    p_dock.add_argument("--library", help="Path to ligand library (.sdf)")
    add_common_wf_args(p_dock)

    # 4. rescore
    p_rescore = subparsers.add_parser("rescore", help="Run MM-GBSA and CNN rescoring")
    add_common_wf_args(p_rescore)

    # 5. md-setup
    p_mds = subparsers.add_parser("md-setup", help="Build MD systems for chosen hits")
    p_mds.add_argument("--hits", default="native", help="Comma-separated hit titles (e.g. 'native,ZINC01')")
    add_common_wf_args(p_mds)

    # 6. md-run
    p_mdr = subparsers.add_parser("md-run", help="Execute chunked MD simulations on Slurm GPUs")
    p_mdr.add_argument("--hits", help="Comma-separated hits to simulate")
    add_common_wf_args(p_mdr)

    # 7. md-analyze
    p_mda = subparsers.add_parser("md-analyze", help="Run trajectory PBC fix, metrics, and MM-PBSA")
    p_mda.add_argument("--traj", help="Direct path to custom .xtc file")
    p_mda.add_argument("--topol", help="Direct path to custom topol.top or PDB")
    add_common_wf_args(p_mda)

    # 8. all
    p_all = subparsers.add_parser("all", help="Execute end-to-end pipeline")
    p_all.add_argument("--hits", help="Comma-separated hit titles for MD")
    add_common_wf_args(p_all)

    # 9. clean
    p_clean = subparsers.add_parser("clean", help="Clean /shared/work scratch directory for current project")

    # 10. status
    p_stat = subparsers.add_parser("status", help="Inspect Slurm queue and Nextflow job status")

    args = parser.parse_args()

    if not args.subcommand:
        parser.print_help()
        sys.exit(0)

    if args.subcommand == "init":
        cmd_init(args)
    elif args.subcommand == "status":
        cmd_status(args)
    elif args.subcommand == "clean":
        cmd_clean(args)
    elif args.subcommand == "prep":
        sys.exit(run_workflow("prepare", args))
    elif args.subcommand == "dock":
        sys.exit(run_workflow("screen", args))
    elif args.subcommand == "md-setup" or args.subcommand == "md-run":
        sys.exit(run_workflow("dynamics", args))
    elif args.subcommand == "md-analyze":
        sys.exit(run_workflow("dynamics", args))
    elif args.subcommand == "all":
        sys.exit(run_workflow("all", args))

if __name__ == "__main__":
    main()

