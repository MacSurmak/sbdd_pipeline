"""
File: cli/commands/traj.py
Description: CLI handler for MD trajectory clustering. Dispatches to Nextflow and auto-registers resulting frames in DB.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Receptor, Trajectory
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_cluster(args):
    cwd = Path.cwd()
    session = init_db()

    traj_name = args.traj
    traj = session.query(Trajectory).filter_by(name=traj_name).first()

    if not traj:
        console.print(
            f"[bold red]FATAL: Trajectory '{traj_name}' not found in DB.[/bold red]"
        )
        sys.exit(1)

    parent_rec = session.query(Receptor).filter_by(name=traj.receptor_name).first()
    if not parent_rec:
        console.print(
            f"[bold red]FATAL: Parent receptor '{traj.receptor_name}' not found.[/bold red]"
        )
        sys.exit(1)

    # Build Plan
    plan_file = cwd / "cluster_plan.json"
    plan = {
        "traj_name": traj.name,
        "xtc": traj.path_xtc,
        "tpr": traj.path_tpr,
        "ref_pdb": parent_rec.path,
        "cutoff": args.cutoff,
        "method": args.method,
    }
    with open(plan_file, "w") as f:
        json.dump(plan, f, indent=2)

    console.print(
        f"[bold cyan]Clustering trajectory '{traj.name}' (Cutoff: {args.cutoff} nm)...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--cluster_plan", str(plan_file.resolve())])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow clustering failed.[/bold red]")
        sys.exit(ret_code)

    # Parse Receipt
    receipt_file = cwd / "cluster_receipt.tsv"
    if not receipt_file.exists():
        console.print("[bold red]FATAL: cluster_receipt.tsv not found![/bold red]")
        sys.exit(1)

    lines = receipt_file.read_text().strip().splitlines()
    count = 0
    for line in lines:
        if not line.strip():
            continue
        cluster_id, pdb_path = line.split("\t")

        # Publish path
        pub_path = cwd / "prepared" / Path(pdb_path).name
        if not pub_path.exists():
            pub_path = Path(pdb_path)

        ex = session.query(Receptor).filter_by(name=cluster_id).first()
        if ex:
            ex.path = str(pub_path)
        else:
            # Auto register as LIGHT prep mode!
            new_rec = Receptor(
                name=cluster_id,
                path=str(pub_path),
                status="unprepared",
                prep_mode="light",
                parent_name=parent_rec.name,
            )
            # Inherit native ligand for docking
            new_rec.meta_dict = {
                "native_ligand": parent_rec.meta_dict.get("native_ligand", "")
            }
            session.add(new_rec)
        count += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Clustering complete! Registered {count} frames as Receptors (light mode).[/bold green]"
    )
    console.print(
        "[dim]Use 'sbdd tree' to view them, and 'sbdd dock' to cross-dock into them.[/dim]"
    )
