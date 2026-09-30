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

    traj_names = [t.strip() for t in args.traj.split(",") if t.strip()]
    trajs = session.query(Trajectory).filter(Trajectory.name.in_(traj_names)).all()

    if len(trajs) != len(traj_names):
        found = {t.name for t in trajs}
        missing = set(traj_names) - found
        console.print(
            f"[bold red]FATAL: Trajectories not found in DB: {missing}[/bold red]"
        )
        sys.exit(1)

    # Build Plan for all requested trajectories
    plan_list = []
    parent_map = {}
    for traj in trajs:
        parent_rec = session.query(Receptor).filter_by(name=traj.receptor_name).first()
        if not parent_rec:
            console.print(
                f"[bold red]FATAL: Parent receptor '{traj.receptor_name}' not found.[/bold red]"
            )
            sys.exit(1)
        parent_map[traj.name] = parent_rec

        plan_list.append(
            {
                "traj_name": traj.name,
                "xtc": traj.path_xtc,
                "tpr": traj.path_tpr,
                "cutoff": args.cutoff,
                "method": args.method,
            }
        )

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Clustering {len(trajs)} trajectory(ies) in parallel on cluster (Cutoff: {args.cutoff} nm)...[/bold cyan]"
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

        # Derive parent traj name from prefix (e.g. 7UA1_300ns_c1 -> 7UA1_300ns)
        traj_prefix = "_".join(cluster_id.split("_")[:-1])
        parent_rec = parent_map.get(traj_prefix)
        parent_name = parent_rec.name if parent_rec else None
        native_lig = parent_rec.meta_dict.get("native_ligand", "") if parent_rec else ""

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
                parent_name=parent_name,
            )
            # Inherit native ligand for docking
            new_rec.meta_dict = {"native_ligand": native_lig}
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
