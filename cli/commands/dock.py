"""
File: cli/commands/dock.py
Description: CLI handler for 'sbdd dock'. Reads DB to construct Cross-Docking plan.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Receptor, Ligand, PoseSet
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_dock(args):
    cwd = Path.cwd()
    session = init_db()

    # 1. Fetch requested Receptors
    rec_names = [n.strip() for n in args.receptors.split(",")]
    receptors = session.query(Receptor).filter(Receptor.name.in_(rec_names)).all()

    if len(receptors) != len(rec_names):
        found = {r.name for r in receptors}
        missing = set(rec_names) - found
        console.print(
            f"[bold red]FATAL: Receptors not found in DB: {missing}[/bold red]"
        )
        sys.exit(1)

    for r in receptors:
        if r.status != "ready":
            console.print(
                f"[bold red]FATAL: Receptor '{r.name}' is not ready! Run 'sbdd prep' first.[/bold red]"
            )
            sys.exit(1)

    # 2. Fetch requested Ligands
    lig_names = [n.strip() for n in args.ligands.split(",")]
    ligands = session.query(Ligand).filter(Ligand.name.in_(lig_names)).all()

    if len(ligands) != len(lig_names):
        found = {l.name for l in ligands}
        missing = set(lig_names) - found
        console.print(f"[bold red]FATAL: Ligands not found in DB: {missing}[/bold red]")
        sys.exit(1)

    # 3. Build Docking Plan
    plan = {"receptors": [], "ligands": []}

    for r in receptors:
        ref_pose_path = "NO_FILE"
        sys_obj = session.query(System).filter_by(receptor_name=r.name).first()
        if sys_obj and sys_obj.pose_name:
            pose_obj = session.query(Pose).filter_by(name=sys_obj.pose_name).first()
            if pose_obj and Path(pose_obj.path).exists():
                ref_pose_path = pose_obj.path

        plan["receptors"].append(
            {
                "name": r.name,
                "path": r.path,
                "ref_pose_path": ref_pose_path,
                "native_ligand": r.meta_dict.get("native_ligand", ""),
            }
        )

    for l in ligands:
        # Determine if this ligand originates from a crystal structure in any system
        has_crystal_pose = (
            session.query(Pose)
            .filter_by(ligand_name=l.name, origin="crystal")
            .first()
            is not None
        )
        plan["ligands"].append(
            {
                "name": l.name,
                "path": l.path,
                "status": l.status,
                "is_native": has_crystal_pose,
            }
        )

    plan_file = cwd / "dock_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan, f, indent=2)

    console.print(
        f"[bold cyan]Cross-Docking Plan generated: {len(receptors)} Receptors × {len(ligands)} Ligand Sets[/bold cyan]"
    )

    # 4. Build Nextflow args
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="screen")
    nf_args.extend(["--dock_plan", str(plan_file.resolve())])

    if getattr(args, "engine", None):
        if "--dock_engine" in nf_args:
            idx = nf_args.index("--dock_engine")
            nf_args[idx + 1] = args.engine
        else:
            nf_args.extend(["--dock_engine", args.engine])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    # 5. Run execution
    ret_code = execute_nextflow(
        "screen", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow docking failed.[/bold red]")
        sys.exit(ret_code)

    # 6. Register results in Database
    poses_name = (
        args.name
        if getattr(args, "name", None)
        else f"screen_{len(receptors)}x{len(ligands)}"
    )
    res_csv = cwd / "reports" / "ranking.csv"
    res_sdf = (
        cwd / "docking" / f"{poses_name}_poses.sdf"
    )  # We will implement SDF merging later if needed, tracking CSV for now

    if res_csv.exists():
        ex_pose = session.query(PoseSet).filter_by(name=poses_name).first()
        if ex_pose:
            ex_pose.path_csv = str(res_csv)
        else:
            new_pose = PoseSet(
                name=poses_name,
                path_sdf="multiple_files",
                path_csv=str(res_csv),
                status="ready",
            )
            session.add(new_pose)
        session.commit()
        console.print(
            f"[bold green]✓ Docking complete! Results registered as PoseSet '{poses_name}'.[/bold green]"
        )

    session.close()
