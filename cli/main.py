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
    """Initializes editable MD protocol playlist and MDP files from assets/mdp templates."""
    cwd = Path.cwd()
    configs = load_project_configs(cwd)
    md_dir = cwd / "configs" / "md_protocol"
    md_dir.mkdir(parents=True, exist_ok=True)

    is_mem = "true" if configs.get("membrane", {}).get("enabled", False) else "false"
    tpl_dir = FRAMEWORK_DIR / "assets" / "mdp"
    script = FRAMEWORK_DIR / "bin" / "generate_mdp.py"

    cmd = [
        sys.executable,
        str(script),
        "--templates_dir",
        str(tpl_dir),
        "--outdir",
        str(md_dir),
        "--is_membrane",
        is_mem,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        console.print(f"[bold red]Failed to generate MD protocol: {res.stderr}[/bold red]")
        sys.exit(1)

    console.print(f"[bold green]✓ MD Protocol initialized from assets/mdp into {md_dir}[/bold green]")
    console.print(
        "[dim]You can now edit the generated .mdp files or reorder protocol.txt before running MD.[/dim]"
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

    # 0. add (Database Registration)
    p_add = subparsers.add_parser(
        "add", help="Register external files into the project database"
    )
    add_subs = p_add.add_subparsers(dest="add_type", required=True)

    # 0.1 add structure / receptor
    p_add_rec = add_subs.add_parser(
        "structure", aliases=["receptor", "complex"], help="Import macromolecular structure(s) into a Target"
    )
    p_add_rec.add_argument(
        "path", nargs="+", help="Path(s) or glob to structure file(s) (.pdb, .cif, .mae)"
    )
    p_add_rec.add_argument(
        "--target", help="Parent target name (e.g. RyR1, TSPO). Defaults to file stem."
    )
    p_add_rec.add_argument(
        "--name", help="Unique identifier for this structural frame"
    )
    p_add_rec.add_argument(
        "--cofactors", nargs="*", default=[], help="Cofactor residue names to retain (e.g. ATP MG)"
    )
    p_add_rec.add_argument(
        "--ligand", help="Native ligand residue name inside structure (e.g. KVR, BTN)"
    )

    # 0.2 add ligand
    p_add_lig = add_subs.add_parser(
        "ligand", help="Register a ligand library or individual molecule"
    )
    p_add_lig.add_argument("path", help="Path to the ligand file (.sdf)")
    p_add_lig.add_argument(
        "--name", help="Identifier for ligand (optional for multi-molecule SDF)"
    )
    p_add_lig.add_argument(
        "--force", action="store_true", help="Overwrite if name exists"
    )

    # 0.3 add traj
    p_add_traj = add_subs.add_parser(
        "traj", help="Register an existing MD trajectory replicate linked to a System"
    )
    p_add_traj.add_argument(
        "--system",
        required=True,
        help="Name of the System in DB (e.g. 9ZKE_crystal_holo)",
    )
    p_add_traj.add_argument(
        "--replicate",
        type=int,
        default=1,
        help="Replicate index (1, 2, 3...)",
    )
    p_add_traj.add_argument(
        "--raw-dir",
        help="Path to directory containing raw chunk files md_*.xtc and sys_ref.tpr",
    )
    p_add_traj.add_argument(
        "--xtc",
        help="Path to raw .xtc file or primary chunk (optional if --raw-dir is provided)",
    )
    p_add_traj.add_argument(
        "--tpr",
        required=True,
        help="Path to .tpr file (e.g. sys_ref.tpr)",
    )
    p_add_traj.add_argument(
        "--length_ns",
        type=float,
        default=0.0,
        help="Total simulation length in nanoseconds (e.g. 100.0)",
    )
    p_add_traj.add_argument(
        "--name",
        help="Unique name for this trajectory (default: <system>_rep<replicate>)",
    )
    p_add_traj.add_argument(
        "--clean-xtc",
        help="Explicit path to pre-cleaned .xtc (optional, if already cleaned)",
    )
    p_add_traj.add_argument(
        "--ref-pdb",
        help="Explicit path to dry reference .pdb (optional, if already cleaned)",
    )
    p_add_traj.add_argument("--force", action="store_true", help="Overwrite if exists")
    # 0.4 tree
    p_tree = subparsers.add_parser(
        "tree", aliases=["ls"], help="Show database lineage tree"
    )

    # 0.4 rm
    p_rm = subparsers.add_parser("rm", help="Remove an item from the database by name")
    p_rm.add_argument("target", help="Name of the item to remove")

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
        "prep",
        aliases=["prepare"],
        help="Prepare raw receptor frames and in-situ poses into ready systems",
    )
    p_prep.add_argument(
        "--target",
        help="Prepare all raw frames for a specific Target",
    )
    p_prep.add_argument(
        "--frames",
        help="Comma-separated list of specific ReceptorFrame names to prepare",
    )
    add_common_wf_args(p_prep)

    # 3. dock / screen
    p_dock = subparsers.add_parser(
        "dock", aliases=["screen"], help="Run distributed docking into ReceptorFrames"
    )
    p_dock.add_argument(
        "--target", help="Parent target name to dock into its ready frames"
    )
    p_dock.add_argument(
        "--frames", "--receptors", dest="frames",
        help="Comma-separated list of specific ready ReceptorFrame names",
    )
    p_dock.add_argument(
        "--ligands",
        help="Comma-separated list of ligand names from DB (optional in --self / --redock mode)",
    )
    p_dock.add_argument(
        "--self", "--redock",
        dest="self_dock",
        action="store_true",
        help="Run 1-to-1 self-docking (redocking) of each receptor with its own native ligand without cross-docking matrix",
    )
    p_dock.add_argument("--name", help="Name for the output PoseSet in DB")
    p_dock.add_argument(
        "--engine",
        choices=["gnina", "glide", "both"],
        default=None,
        help="Docking engine choice (default: from configs)",
    )
    add_common_wf_args(p_dock)

    # 3.1 select-hits
    p_hits = subparsers.add_parser(
        "select-hits",
        help="Select top candidate hits from docking rankings and promote to Systems",
    )
    p_hits.add_argument(
        "--poses", help="Name of PoseSet in DB to select hits from (default: latest)"
    )
    p_hits.add_argument(
        "--top", type=int, default=5, help="Number of top hits to extract (default: 5)"
    )
    p_hits.add_argument(
        "--system-prefix",
        default="hit",
        help="Prefix for generated Holo systems (default: hit)",
    )

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
        help="Execute chunked MD simulations for a System using protocol.txt playlist",
    )
    p_dyn.add_argument(
        "--system", help="Name of the physical System to simulate (from 'sbdd tree')"
    )
    p_dyn.add_argument(
        "--replicates",
        type=int,
        help="Number of independent replicates (default: from config)",
    )
    p_dyn.add_argument(
        "--total_ns",
        type=float,
        help="Total simulation duration in ns (default: from config)",
    )
    add_common_wf_args(p_dyn)

    # 7. traj-mmpbsa
    p_mmpbsa = subparsers.add_parser(
        "traj-mmpbsa", help="Run standalone MM-PBSA free energy calculation on clean trajectories"
    )
    p_mmpbsa.add_argument(
        "--systems", help="Limit to specific System names"
    )
    p_mmpbsa.add_argument(
        "--target", help="Limit to specific Target"
    )
    p_mmpbsa.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    add_common_wf_args(p_mmpbsa)

    # 7. md-analyze
    p_mda = subparsers.add_parser(
        "md-analyze", help="Run trajectory PBC fix, metrics, and MM-PBSA"
    )
    p_mda.add_argument("--traj", help="Direct path to custom .xtc file")
    p_mda.add_argument("--tpr", help="Direct path to custom .tpr file")
    p_mda.add_argument("--topol", help="Direct path to custom topol.top or PDB")
    p_mda.add_argument("--ndx", help="Direct path to custom index.ndx file")
    add_common_wf_args(p_mda)

    # 7.0 traj-flag
    p_flag = subparsers.add_parser(
        "traj-flag",
        help="Flag or disable a trajectory (e.g. ligand drifted/unbound) to exclude from clustering",
    )
    p_flag.add_argument("trajectory", help="Name of the Trajectory record in DB")
    p_flag.add_argument(
        "--status",
        choices=["cleaned", "drifted", "failed", "disabled", "raw"],
        default="drifted",
        help="Target QC status (default: drifted)",
    )
    p_flag.add_argument("--disable", action="store_true", help="Shortcut to set status='disabled'")
    p_flag.add_argument("--enable", action="store_true", help="Shortcut to restore status='cleaned'")

    # 7.1 traj-clean
    p_clean_tr = subparsers.add_parser(
        "traj-clean",
        help="Remove solvent, fix PBC, center on protein, and fit trajectory for all uncleaned runs",
    )
    p_clean_tr.add_argument(
        "--systems",
        help="Limit to specific System names",
    )
    p_clean_tr.add_argument(
        "--target",
        help="Limit to specific Target",
    )
    p_clean_tr.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    add_common_wf_args(p_clean_tr)

    # 7.2 traj-metrics
    p_metrics = subparsers.add_parser(
        "traj-metrics",
        help="Calculate RMSD (protein, ligand, and cofactors), RMSF, and contacts for cleaned trajectories",
    )
    p_metrics.add_argument(
        "--systems", help="Limit to specific System names"
    )
    p_metrics.add_argument(
        "--target", help="Limit to specific Target"
    )
    p_metrics.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    add_common_wf_args(p_metrics)

    # 8. traj-cluster
    p_cluster = subparsers.add_parser(
        "traj-cluster",
        help="Cluster clean MD trajectory ensembles and automatically register frames in DB for docking",
    )
    p_cluster.add_argument(
        "--systems",
        help="Limit clustering to specific System names",
    )
    p_cluster.add_argument(
        "--target",
        help="Limit clustering to a specific Target",
    )
    p_cluster.add_argument(
        "--traj", help="Explicit comma-separated list of Trajectory names"
    )
    p_cluster.add_argument(
        "--cutoff", type=float, default=0.15, help="RMSD cutoff in nm (default: 0.15)"
    )
    p_cluster.add_argument(
        "--dt",
        type=int,
        default=100,
        help="Frame stride in ps to accelerate clustering (default: 100 ps)",
    )
    p_cluster.add_argument(
        "--group",
        default="COMPLEX",
        help="Group for RMSD distance matrix (default: COMPLEX)",
    )
    p_cluster.add_argument(
        "--method", default="gromos", help="Clustering method (default: gromos)"
    )
    p_cluster.add_argument(
        "--max-clusters",
        type=int,
        default=5,
        help="Maximum number of top-populated clusters to retain (default: 5)",
    )
    add_common_wf_args(p_cluster)

    # 9. clean
    p_clean = subparsers.add_parser(
        "clean", help="Clean /shared/work scratch directory for current project"
    )

    # 9. status
    p_stat = subparsers.add_parser(
        "status", help="Inspect Slurm queue and Nextflow job status"
    )

    # 10. watch / progress
    p_watch = subparsers.add_parser(
        "watch",
        aliases=["progress", "top"],
        help="Monitor real-time task progress via cluster UDP telemetry",
    )
    p_watch.add_argument(
        "--port", type=int, default=8089, help="Telemetry UDP port (default: 8089)"
    )

    args = parser.parse_args()

    if not args.subcommand:
        parser.print_help()
        sys.exit(0)

    subcmd = args.subcommand

    # Init DB unconditionally to ensure it exists on any command run inside a project
    if subcmd not in ["init", "help"]:
        from cli.core.db import init_db

        init_db()

    if subcmd == "add":
        from cli.commands.add import add_receptor, add_ligand, add_traj

        if args.add_type in ["structure", "receptor", "complex"]:
            add_receptor(args)
        elif args.add_type == "ligand":
            add_ligand(args)
        elif args.add_type == "traj":
            add_traj(args)
    elif subcmd in ["tree", "ls"]:
        from cli.commands.tree import cmd_tree

        cmd_tree(args)
    elif subcmd == "rm":
        from cli.commands.rm import cmd_rm

        cmd_rm(args)
    elif subcmd == "init":
        cmd_init(args)
    elif subcmd == "md-init":
        cmd_md_init(args)
    elif subcmd == "status":
        cmd_status(args)
    elif subcmd == "clean":
        cmd_clean(args)
    elif subcmd in ["prep", "prepare"]:
        from cli.commands.prep import cmd_prep

        cmd_prep(args)
    elif subcmd in ["dock", "screen"]:
        from cli.commands.dock import cmd_dock

        cmd_dock(args)
    elif subcmd == "select-hits":
        from cli.commands.select_hits import cmd_select_hits

        cmd_select_hits(args)
    elif subcmd in ["dynamics", "md", "md-run", "md-setup"]:
        from cli.commands.dynamics import cmd_dynamics

        cmd_dynamics(args)
    elif subcmd == "md-analyze":
        sys.exit(run_workflow("dynamics", args))
    elif subcmd == "traj-flag":
        from cli.commands.traj import cmd_flag_traj

        cmd_flag_traj(args)
    elif subcmd == "traj-clean":
        from cli.commands.traj import cmd_clean_traj

        cmd_clean_traj(args)
    elif subcmd == "traj-mmpbsa":
        from cli.commands.traj import cmd_mmpbsa

        cmd_mmpbsa(args)
    elif subcmd == "traj-metrics":
        from cli.commands.traj import cmd_metrics

        cmd_metrics(args)
    elif subcmd == "traj-cluster":
        from cli.commands.traj import cmd_cluster

        cmd_cluster(args)
    elif subcmd in ["watch", "progress", "top"]:
        from cli.core.monitor import run_monitor

        run_monitor(port=args.port)


if __name__ == "__main__":
    main()
