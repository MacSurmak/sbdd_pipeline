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
    console.print(
        "  1. Copy your structures (.cif or .pdb) into the project directory."
    )
    console.print("  2. Edit configs/manifest.yaml to describe your targets.")
    console.print("  3. Review configurations (e.g., md.yaml, gnina.yaml).")
    console.print(
        "  4. Run [bold cyan]sbdd prep[/bold cyan] to standardize the structures!\n"
    )


def cmd_status(args):
    """Inspects active Slurm jobs and Nextflow status."""
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

    # Check Nextflow log in cwd
    log_f = Path.cwd() / ".nextflow.log"
    if log_f.exists():
        console.rule("[bold cyan]Recent Nextflow Activity[/bold cyan]")
        tail = subprocess.run(
            ["tail", "-n", "15", str(log_f)], capture_output=True, text=True
        )
        console.print(tail.stdout)


def cmd_md_init(args):
    """Initializes the editable MD protocol playlist and MDP files."""
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    md_dir = cwd / "configs" / "md_protocol"
    md_dir.mkdir(parents=True, exist_ok=True)

    is_mem = configs.get("membrane", {}).get("enabled", False)

    # Generate clean MDP templates with decoupled protein and ligand preprocessor controls
    content_minim = "define = -DPOSRES -DPOSRES_LIG\nintegrator = steep\nemtol = 500.0\nemstep = 0.01\nnsteps = 50000\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nvdwtype = cutoff\nrvdw = 1.0\npbc = xyz\nconstraints = h-bonds\nconstraint-algorithm = lincs\n"
    content_nvt = "define = -DPOSRES -DPOSRES_LIG\nintegrator = md\ndt = 0.001\nnsteps = 50000\nnstxout-compressed = 5000\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nvdwtype = cutoff\nrvdw = 1.0\npbc = xyz\ntcoupl = V-rescale\ntc-grps = PROT_LIG SOL_ION\ntau-t = 0.1 0.1\nref-t = 310.15 310.15\npcoupl = no\nconstraints = h-bonds\nconstraint-algorithm = lincs\ncontinuation = no\ngen_vel = yes\ngen_temp = 310.15\n"
    content_npt = "define = -DPOSRES_250 -DPOSRES_LIG_250\nintegrator = md\ndt = 0.002\nnsteps = 100000\nnstxout-compressed = 5000\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nvdwtype = cutoff\nrvdw = 1.0\npbc = xyz\ntcoupl = V-rescale\ntc-grps = PROT_LIG SOL_ION\ntau-t = 1.0 1.0\nref-t = 310.15 310.15\npcoupl = Parrinello-Rahman\npcoupltype = isotropic\ntau-p = 5.0\nref-p = 1.0\ncompressibility = 4.5e-5\nrefcoord-scaling = com\nconstraints = h-bonds\ncontinuation = yes\n"
    content_md = "integrator = md\ndt = 0.002\nnsteps = 500000\nnstxout-compressed = 5000\ncutoff-scheme = Verlet\ncoulombtype = PME\nrcoulomb = 1.0\nvdwtype = cutoff\nrvdw = 1.0\npbc = xyz\ntcoupl = V-rescale\ntc-grps = PROT_LIG SOL_ION\ntau-t = 1.0 1.0\nref-t = 310.15 310.15\npcoupl = Parrinello-Rahman\npcoupltype = isotropic\ntau-p = 5.0\nref-p = 1.0\ncompressibility = 4.5e-5\nconstraints = h-bonds\ncontinuation = yes\n"

    if is_mem:
        content_npt = (
            content_npt.replace("isotropic", "semiisotropic")
            .replace("ref-p = 1.0", "ref-p = 1.0 1.0")
            .replace("compressibility = 4.5e-5", "compressibility = 4.5e-5 4.5e-5")
        )
        content_md = (
            content_md.replace("isotropic", "semiisotropic")
            .replace("ref-p = 1.0", "ref-p = 1.0 1.0")
            .replace("compressibility = 4.5e-5", "compressibility = 4.5e-5 4.5e-5")
        )

    (md_dir / "minim.mdp").write_text(content_minim)
    (md_dir / "nvt.mdp").write_text(content_nvt)
    (md_dir / "npt.mdp").write_text(content_npt)
    (md_dir / "md.mdp").write_text(content_md)

    (md_dir / "protocol.txt").write_text("minim.mdp\nnvt.mdp\nnpt.mdp\n")

    console.print(f"[bold green]✓ MD Protocol initialized in {md_dir}[/bold green]")
    console.print(
        "[dim]You can now manually edit the .mdp files and reorder protocol.txt before running MD.[/dim]"
    )


def cmd_clean(args):
    """Cleans the project's scratch directory."""
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    proj_name = configs.get("project", {}).get("name", cwd.name)
    workdir = configs.get("project", {}).get("workdir", f"/shared/work/{proj_name}")
    clean_shared_scratch(workdir, verbose=True)


def run_workflow(entrypoint: str, args):
    """Generic dispatcher from CLI flags and configs to Nextflow execution."""
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint=entrypoint)

    # CLI flag overrides
    if getattr(args, "receptor", None):
        nf_args.extend(["--receptor", str(Path(args.receptor).resolve())])
    if getattr(args, "library", None):
        nf_args.extend(["--library", str(Path(args.library).resolve())])
    if getattr(args, "engine", None):
        nf_args.extend(["--dock_engine", args.engine])
    if getattr(args, "hits", None):
        nf_args.extend(["--md_hits", args.hits])
    if getattr(args, "traj", None):
        nf_args.extend(["--traj", str(Path(args.traj).resolve())])
    if getattr(args, "tpr", None):
        nf_args.extend(["--tpr", str(Path(args.tpr).resolve())])
    if getattr(args, "topol", None):
        nf_args.extend(["--topol", str(Path(args.topol).resolve())])
    if getattr(args, "ndx", None):
        nf_args.extend(["--ndx", str(Path(args.ndx).resolve())])
    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)
    auto_clean = getattr(args, "clean_scratch", False)

    return execute_nextflow(
        entrypoint,
        nf_args,
        cwd,
        profile=profile,
        resume=resume,
        auto_clean=auto_clean,
        configs=configs,
    )


def main():
    parser = argparse.ArgumentParser(
        prog="sbdd", description="SBDD Pipeline Nextflow HPC Orchestrator CLI"
    )
    subparsers = parser.add_subparsers(dest="subcommand", help="Available subcommands")

    # 1. init
    p_init = subparsers.add_parser("init", help="Initialize new project workspace")
    p_init.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Target project directory (default: current)",
    )
    p_init.add_argument(
        "--force", action="store_true", help="Overwrite existing templates"
    )

    # Common options for workflows
    def add_common_wf_args(p):
        p.add_argument("--local", action="store_true", help="Run locally without Slurm")
        p.add_argument(
            "--no-resume",
            action="store_true",
            help="Start from scratch without -resume",
        )
        p.add_argument(
            "--clean-scratch",
            action="store_true",
            help="Delete /shared/work upon completion (default: False)",
        )
        p.add_argument(
            "--urgent",
            action="store_true",
            help="Route tasks to urgent Slurm partition",
        )

    # 2. prep / prepare
    p_prep = subparsers.add_parser(
        "prep", aliases=["prepare"], help="Run receptor & ligand preparation"
    )
    p_prep.add_argument("--receptor", help="Path to input structure (.cif/.pdb)")
    p_prep.add_argument("--library", help="Path to ligand library (.sdf)")
    add_common_wf_args(p_prep)

    # 3. dock / screen
    p_dock = subparsers.add_parser(
        "dock", aliases=["screen"], help="Run distributed docking across GPU nodes"
    )
    p_dock.add_argument(
        "--engine",
        choices=["gnina", "glide", "both"],
        default=None,
        help="Docking engine choice (default: from configs)",
    )
    p_dock.add_argument("--library", help="Path to ligand library (.sdf)")
    add_common_wf_args(p_dock)

    # 4. rescore
    p_rescore = subparsers.add_parser("rescore", help="Run MM-GBSA and CNN rescoring")
    add_common_wf_args(p_rescore)

    # 5. md-init
    p_mdinit = subparsers.add_parser(
        "md-init", help="Initialize editable MD protocol playlist"
    )

    # 6. md / dynamics
    p_dyn = subparsers.add_parser(
        "dynamics",
        aliases=["md", "md-run", "md-setup"],
        help="Execute chunked MD simulations using protocol.txt playlist",
    )
    add_common_wf_args(p_dyn)

    # 7. md-analyze
    p_mda = subparsers.add_parser(
        "md-analyze", help="Run trajectory PBC fix, metrics, and MM-PBSA"
    )
    p_mda.add_argument("--traj", help="Direct path to custom .xtc file")
    p_mda.add_argument("--tpr", help="Direct path to custom .tpr file")
    p_mda.add_argument("--topol", help="Direct path to custom topol.top or PDB")
    p_mda.add_argument("--ndx", help="Direct path to custom index.ndx file")
    add_common_wf_args(p_mda)

    # 7. clean
    p_clean = subparsers.add_parser(
        "clean", help="Clean /shared/work scratch directory for current project"
    )

    # 9. status
    p_stat = subparsers.add_parser(
        "status", help="Inspect Slurm queue and Nextflow job status"
    )

    args = parser.parse_args()

    if not args.subcommand:
        parser.print_help()
        sys.exit(0)

    subcmd = args.subcommand
    if subcmd == "init":
        cmd_init(args)
    elif subcmd == "md-init":
        cmd_md_init(args)
    elif subcmd == "status":
        cmd_status(args)
    elif subcmd == "clean":
        cmd_clean(args)
    elif subcmd in ["prep", "prepare"]:
        sys.exit(run_workflow("prepare", args))
    elif subcmd in ["dock", "screen"]:
        sys.exit(run_workflow("screen", args))
    elif subcmd in ["dynamics", "md", "md-run", "md-setup", "md-analyze"]:
        sys.exit(run_workflow("dynamics", args))


if __name__ == "__main__":
    main()
