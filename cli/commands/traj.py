"""
File: cli/commands/traj.py
Description: CLI handler for MD trajectory clustering. Dispatches to Nextflow and auto-registers resulting frames in DB.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Receptor, Trajectory, System, Pose
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def _query_trajectories(session, args):
    """Helper to query trajectories linked to Systems or Receptors."""
    if getattr(args, "all", False):
        trajs = session.query(Trajectory).all()
        if not trajs:
            console.print("[bold red]FATAL: No trajectories registered in database.[/bold red]")
            sys.exit(1)
        return trajs
    elif getattr(args, "systems", None):
        target_sys = [s.strip() for s in args.systems.split(",") if s.strip()]
        trajs = session.query(Trajectory).filter(Trajectory.system_name.in_(target_sys)).all()
        if not trajs:
            console.print(f"[bold red]FATAL: No trajectories found for system(s): {target_sys}[/bold red]")
            sys.exit(1)
        return trajs
    elif getattr(args, "receptors", None):
        target_recs = [r.strip() for r in args.receptors.split(",") if r.strip()]
        systems = session.query(System).filter(System.receptor_name.in_(target_recs)).all()
        sys_names = [s.name for s in systems]
        trajs = session.query(Trajectory).filter(Trajectory.system_name.in_(sys_names)).all()
        if not trajs:
            console.print(f"[bold red]FATAL: No trajectories found for receptor(s): {target_recs}[/bold red]")
            sys.exit(1)
        return trajs
    elif getattr(args, "traj", None):
        traj_names = [t.strip() for t in args.traj.split(",") if t.strip()]
        trajs = session.query(Trajectory).filter(Trajectory.name.in_(traj_names)).all()
        if len(trajs) != len(traj_names):
            found = {t.name for t in trajs}
            missing = set(traj_names) - found
            console.print(f"[bold red]FATAL: Trajectories not found in DB: {missing}[/bold red]")
            sys.exit(1)
        return trajs
    else:
        console.print("[bold red]FATAL: Please specify --systems, --receptors, --all, or --traj[/bold red]")
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

        sys_obj = session.query(System).filter_by(name=t.system_name).first()
        rec_name = sys_obj.receptor_name if sys_obj else t.system_name
        rec_id = rec_name.replace("_ready", "")
        xtc_path = Path(t.path_xtc).resolve()

        # Destination inside project directory next to replicate
        if xtc_path.parent.name == "raw":
            clean_dir = xtc_path.parent.parent / "clean"
        else:
            clean_dir = cwd / "md" / rec_id / f"rep_{rep_idx}" / "clean"

        plan_list.append(
            {
                "name": t.name,
                "rec_id": rec_id,
                "rep": rep_idx,
                "xtc": str(xtc_path),
                "tpr": str(Path(t.path_tpr).resolve()),
                "clean_dir": str(clean_dir.resolve()),
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

    for item in plan_list:
        t_record = session.query(Trajectory).filter_by(name=item["name"]).first()
        if t_record:
            clean_dir = Path(item["clean_dir"])
            t_record.path_clean_xtc = str(clean_dir / "trajectory_clean.xtc")
            t_record.path_ref_pdb = str(clean_dir / "complex_ref.pdb")
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
    from rich.table import Table

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

        sys_obj = session.query(System).filter_by(name=t.system_name).first()
        parent_rec = (
            session.query(Receptor).filter_by(name=sys_obj.receptor_name).first()
            if sys_obj
            else None
        )
        native_lig = parent_rec.meta_dict.get("native_ligand", "") if parent_rec else ""
        cofactors = parent_rec.meta_dict.get("cofactors", []) if parent_rec else []
        rec_name = sys_obj.receptor_name if sys_obj else t.system_name
        rec_id = rec_name.replace("_ready", "")

        # Route metrics directly into replicate folder
        clean_xtc_path = Path(t.path_clean_xtc).resolve()
        if clean_xtc_path.parent.name == "clean":
            metrics_dir = clean_xtc_path.parent.parent / "metrics"
        else:
            metrics_dir = cwd / "md" / rec_id / f"rep_{rep_idx}" / "metrics"

        plan_list.append(
            {
                "name": t.name,
                "rec_id": rec_id,
                "rep": rep_idx,
                "clean_xtc": t.path_clean_xtc,
                "ref_pdb": t.path_ref_pdb,
                "ligand": native_lig,
                "cofactors": cofactors,
                "metrics_dir": str(metrics_dir.resolve()),
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

    # Render summary table in console
    table = Table(title="Molecular Dynamics Ensemble Stability Summary", show_header=True)
    table.add_column("Trajectory", style="bold cyan")
    table.add_column("Prot BB RMSD (Å)", justify="center")
    table.add_column("Ligand RMSD (Å)", justify="center", style="green")
    table.add_column("Cofactors RMSD (Å)", justify="center", style="yellow")
    table.add_column("Plot Deliverable", style="dim")

    for item in plan_list:
        summary_f = Path(item["metrics_dir"]) / "md_metrics_summary.json"
        plot_f = Path(item["metrics_dir"]) / "rmsd_plot.png"
        if summary_f.exists():
            data = json.loads(summary_f.read_text())
            p_mean = data.get("Prot_BB_RMSD_Mean_A", "-")
            p_std = data.get("Prot_BB_RMSD_Std_A", "-")
            prot_str = f"{p_mean} ± {p_std}" if p_mean != "-" else "-"

            lig_str = "-"
            cof_strs = []
            mols = data.get("Molecules", {})
            for m_key, m_val in mols.items():
                if m_key.startswith("Lig_"):
                    lig_str = f"{m_val.get('Mean_A', '-')} ± {m_val.get('Std_A', '-')}"
                elif m_key.startswith("Cofactor_"):
                    c_name = m_key.replace("Cofactor_", "").replace("_RMSD_A", "")
                    cof_strs.append(f"{c_name}: {m_val.get('Mean_A', '-')}")

            cof_display = ", ".join(cof_strs) if cof_strs else "None"
            rel_plot = str(plot_f.relative_to(cwd)) if plot_f.exists() else "-"
            table.add_row(item["name"], prot_str, lig_str, cof_display, rel_plot)

    console.print(table)
    session.close()


def cmd_cluster(args):
    """CLI handler for 'sbdd traj-cluster': clusters clean trajectory ensembles across Systems."""
    cwd = Path.cwd()
    session = init_db()
    trajs = _query_trajectories(session, args)

    # 1. Verify all trajectories are cleaned
    for t in trajs:
        if t.status != "cleaned" or not t.path_clean_xtc or not t.path_ref_pdb:
            console.print(f"[bold red]FATAL: Trajectory '{t.name}' must be cleaned first. Run 'sbdd traj-clean'![/bold red]")
            sys.exit(1)

    # 2. Group requested trajectories by System
    groups = {}
    for t in trajs:
        groups.setdefault(t.system_name, []).append(t)

    plan_list = []
    parent_map = {}

    for sys_name, t_list in groups.items():
        sys_obj = session.query(System).filter_by(name=sys_name).first()
        if not sys_obj:
            console.print(f"[bold red]FATAL: System '{sys_name}' not found.[/bold red]")
            sys.exit(1)

        task_id = f"{sys_name}_ens"
        parent_map[task_id] = sys_obj

        clean_xtcs = [t.path_clean_xtc for t in t_list]
        ref_pdb = t_list[0].path_ref_pdb

        plan_list.append(
            {
                "task_id": task_id,
                "system": sys_name,
                "clean_xtcs": clean_xtcs,
                "ref_pdb": ref_pdb,
                "cutoff": args.cutoff,
                "method": args.method,
                "dt": args.dt,
                "group": args.group,
            }
        )

    plan_file = cwd / "cluster_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Clustering clean ensembles across {len(groups)} System(s) (Cutoff: {args.cutoff} nm, Group: {args.group}, dt: {args.dt} ps)...[/bold cyan]"
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
        parent_sys = parent_map.get(traj_prefix)
        rec_name = parent_sys.receptor_name if parent_sys else None

        ex = session.query(Receptor).filter_by(name=cluster_id).first()
        if ex:
            ex.path = str(pub_path)
        else:
            # Register centroid as a Receptor in light mode
            new_rec = Receptor(
                name=cluster_id,
                path=str(pub_path),
                status="unprepared",
                prep_mode="light",
                parent_name=rec_name,
            )
            session.add(new_rec)

            # Auto-register corresponding cluster System
            cluster_sys_name = f"{cluster_id}_sys"
            new_sys = System(
                name=cluster_sys_name,
                receptor_name=cluster_id,
                pose_name=parent_sys.pose_name if parent_sys else None,
                status="unprepared",
            )
            session.add(new_sys)
        count += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Clustering complete! Registered {count} clean complex frames as Receptors (light mode).[/bold green]"
    )
    console.print(
        "[dim]Protein, Ligand, and Cofactors are preserved in each cluster centroid.[/dim]"
    )
