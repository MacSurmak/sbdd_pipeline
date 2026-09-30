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


def _query_trajectories(session, args):
    """Helper to query trajectories based on CLI flags."""
    if getattr(args, "all", False):
        trajs = session.query(Trajectory).all()
        if not trajs:
            console.print(
                "[bold red]FATAL: No trajectories registered in database.[/bold red]"
            )
            sys.exit(1)
        return trajs
    elif getattr(args, "receptors", None):
        target_recs = [r.strip() for r in args.receptors.split(",") if r.strip()]
        trajs = (
            session.query(Trajectory)
            .filter(Trajectory.receptor_name.in_(target_recs))
            .all()
        )
        if not trajs:
            console.print(
                f"[bold red]FATAL: No trajectories found for receptor(s): {target_recs}[/bold red]"
            )
            sys.exit(1)
        return trajs
    elif getattr(args, "traj", None):
        traj_names = [t.strip() for t in args.traj.split(",") if t.strip()]
        trajs = session.query(Trajectory).filter(Trajectory.name.in_(traj_names)).all()
        if len(trajs) != len(traj_names):
            found = {t.name for t in trajs}
            missing = set(traj_names) - found
            console.print(
                f"[bold red]FATAL: Trajectories not found in DB: {missing}[/bold red]"
            )
            sys.exit(1)
        return trajs
    else:
        console.print(
            "[bold red]FATAL: Please specify --receptors, --all, or --traj[/bold red]"
        )
        sys.exit(1)


def cmd_clean_traj(args):
    """CLI handler for 'sbdd traj-clean': fixes PBC, centers on protein, fits, and extracts COMPLEX."""
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(session, args)

    plan_list = []
    for t in trajs:
        # Resolve index.ndx from parent MD folder if present, or TPR
        tpr_path = Path(t.path_tpr)
        ndx_path = tpr_path.parent / "index.ndx"
        if not ndx_path.exists():
            ndx_path = tpr_path.parent.parent / "index.ndx"
        if not ndx_path.exists():
            ndx_path = tpr_path  # Will be generated dynamically if needed

        # Determine replicate index
        rep_idx = 1
        for part in t.name.split("_"):
            if part.startswith("rep") or part.startswith("R"):
                try:
                    rep_idx = int("".join(filter(str.isdigit, part)))
                except ValueError:
                    pass

        plan_list.append(
            {
                "name": t.name,
                "rep": rep_idx,
                "xtc": t.path_xtc,
                "tpr": t.path_tpr,
            }
        )

    plan_file = cwd / "clean_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Cleaning {len(trajs)} trajectory(ies) (PBC whole -> Center -> Fit -> Extract COMPLEX)...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--clean_plan", str(plan_file.resolve())])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow trajectory cleaning failed.[/bold red]")
        sys.exit(ret_code)

    receipt_file = cwd / "clean_receipt.tsv"
    if not receipt_file.exists():
        console.print("[bold red]FATAL: clean_receipt.tsv not found![/bold red]")
        sys.exit(1)

    for line in receipt_file.read_text().strip().splitlines():
        if not line.strip():
            continue
        tname, clean_xtc, ref_pdb = line.split("\t")
        t_record = session.query(Trajectory).filter_by(name=tname).first()
        if t_record:
            t_record.path_clean_xtc = clean_xtc
            t_record.path_ref_pdb = ref_pdb
            t_record.status = "cleaned"

    session.commit()
    session.close()
    console.print(
        f"\n[bold green]✓ Successfully cleaned {len(trajs)} trajectories and updated the database![/bold green]"
    )
    console.print(
        "[dim]Solvent stripped. Full COMPLEX (Protein + Ligand + Cofactors) retained.[/dim]"
    )


def cmd_metrics(args):
    """CLI handler for 'sbdd traj-metrics': calculates RMSD for protein, ligand, and cofactors."""
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(session, args)

    plan_list = []
    for t in trajs:
        if t.status != "cleaned" or not t.path_clean_xtc or not t.path_ref_pdb:
            console.print(
                f"[bold red]FATAL: Trajectory '{t.name}' is not cleaned yet! Run 'sbdd traj-clean' first.[/bold red]"
            )
            sys.exit(1)

        rep_idx = 1
        for part in t.name.split("_"):
            if part.startswith("rep") or part.startswith("R"):
                try:
                    rep_idx = int("".join(filter(str.isdigit, part)))
                except ValueError:
                    pass

        plan_list.append(
            {
                "name": t.name,
                "rep": rep_idx,
                "clean_xtc": t.path_clean_xtc,
                "ref_pdb": t.path_ref_pdb,
            }
        )

    plan_file = cwd / "metrics_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Calculating metrics for {len(trajs)} clean trajectory(ies) (RMSD Protein + Ligand + Cofactors)...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-analyze")
    nf_args.extend(["--metrics_plan", str(plan_file.resolve())])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "md-analyze", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )
    if ret_code != 0:
        console.print("[bold red]Nextflow metrics calculation failed.[/bold red]")
        sys.exit(ret_code)

    console.print(
        f"\n[bold green]✓ Trajectory metrics calculated successfully![/bold green]"
    )
    session.close()


def cmd_cluster(args):
    """CLI handler for 'sbdd traj-cluster': clusters clean trajectory ensembles (full timeline)."""
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(session, args)

    # 1. Verify all trajectories are cleaned
    for t in trajs:
        if t.status != "cleaned" or not t.path_clean_xtc or not t.path_ref_pdb:
            console.print(
                f"[bold red]FATAL: Trajectory '{t.name}' must be cleaned first. Run 'sbdd traj-clean'![/bold red]"
            )
            sys.exit(1)

    # 2. Group requested trajectories automatically by parent Receptor
    groups = {}
    for t in trajs:
        groups.setdefault(t.receptor_name, []).append(t)

    plan_list = []
    parent_map = {}

    for rec_name, t_list in groups.items():
        parent_rec = session.query(Receptor).filter_by(name=rec_name).first()
        if not parent_rec:
            console.print(
                f"[bold red]FATAL: Parent receptor '{rec_name}' not found.[/bold red]"
            )
            sys.exit(1)

        task_id = f"{rec_name}_ens"
        parent_map[task_id] = parent_rec

        clean_xtcs = [t.path_clean_xtc for t in t_list]
        ref_pdb = t_list[0].path_ref_pdb

        plan_list.append(
            {
                "task_id": task_id,
                "receptor": rec_name,
                "clean_xtcs": clean_xtcs,
                "ref_pdb": ref_pdb,
                "cutoff": args.cutoff,
                "method": args.method,
            }
        )

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Clustering clean ensembles across {len(groups)} Receptor(s) (Cutoff: {args.cutoff} nm, full trajectory)...[/bold cyan]"
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

        pub_path = cwd / "prepared" / Path(pdb_path).name
        if not pub_path.exists():
            pub_path = Path(pdb_path)

        traj_prefix = "_".join(cluster_id.split("_")[:-1])
        parent_rec = parent_map.get(traj_prefix)
        parent_name = parent_rec.name if parent_rec else None
        native_lig = parent_rec.meta_dict.get("native_ligand", "") if parent_rec else ""

        ex = session.query(Receptor).filter_by(name=cluster_id).first()
        if ex:
            ex.path = str(pub_path)
        else:
            new_rec = Receptor(
                name=cluster_id,
                path=str(pub_path),
                status="unprepared",
                prep_mode="light",
                parent_name=parent_name,
            )
            new_rec.meta_dict = {"native_ligand": native_lig}
            session.add(new_rec)
        count += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Clustering complete! Registered {count} clean complex frames as Receptors (light mode).[/bold green]"
    )
    console.print(
        "[dim]Protein, Ligand, and Cofactors are preserved in each cluster centroid.[/dim]"
    )
