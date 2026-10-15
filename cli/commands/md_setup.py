"""
File: cli/commands/md_setup.py
Description: CLI handler for 'sbdd md-setup'.
Gathers 'raw' Systems, generates topologies, and atomically registers them as 'ready'.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Target, System
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def run_md_setup(
    system: str | None = None,
    target: str | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    cwd = Path.cwd()
    session = init_db()
    configs = load_project_configs(cwd)

    query = session.query(System).filter_by(state="raw")
    if system:
        query = query.filter_by(name=system)
    elif target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(f"[bold red]FATAL: Target '{target}' not found.[/bold red]")
            sys.exit(1)
        query = query.filter_by(target_id=target_obj.id)

    raw_systems = query.all()
    if not raw_systems:
        console.print(
            "[green]No raw Systems found to setup. Topologies are up to date![/green]"
        )
        session.close()
        return

    plan_list = []
    for s in raw_systems:
        frame = s.receptor_frame
        if not frame.path_pdb:
            console.print(f"[bold red]FATAL: ReceptorFrame '{frame.name}' missing path_pdb. Run 'sbdd prep' first.[/bold red]")
            sys.exit(1)
        
        rec_pdb = Path(frame.path_pdb).resolve()
        if not rec_pdb.exists():
            raise FileNotFoundError(f"Prepared receptor PDB missing: {rec_pdb}")

        lig_sdf = None
        extra_sdfs = []
        extra_resnames = []
        mol2_paths = []
        frcmod_paths = []

        # 1. Active ligand
        if s.pose:
            if s.pose.substance.substance_type != "ion":
                if not s.pose.path_mol2 or not s.pose.path_frcmod:
                    console.print(f"[bold red]FATAL: System '{s.name}' contains unparametrized ligand '{s.pose.name}'. Run 'sbdd prep'![/bold red]")
                    sys.exit(1)
                mol2_paths.append(str(Path(s.pose.path_mol2).resolve()))
                frcmod_paths.append(str(Path(s.pose.path_frcmod).resolve()))
            
            p_sdf = s.pose.path_sdf or s.pose.path
            if not p_sdf:
                console.print(f"[bold red]FATAL: Pose '{s.pose.name}' missing path_sdf.[/bold red]")
                sys.exit(1)
            lig_sdf = str(Path(p_sdf).resolve())

        # 2. Cofactors and Ions (Separated from pure APO ReceptorFrame)
        for c in s.cofactors:
            c_sdf = c.path_sdf or c.path
            if not c_sdf or not Path(c_sdf).exists():
                console.print(f"[bold red]FATAL: System '{s.name}' cofactor/ion '{c.name}' missing SDF deliverable.[/bold red]")
                sys.exit(1)
            extra_sdfs.append(str(Path(c_sdf).resolve()))
            c_resname = c.res_name
            extra_resnames.append(c_resname)

            if c.substance.substance_type != "ion":
                if not c.path_mol2 or not c.path_frcmod:
                    console.print(f"[bold red]FATAL: System '{s.name}' contains unparametrized cofactor '{c.name}'. Run 'sbdd prep'![/bold red]")
                    sys.exit(1)
                mol2_paths.append(str(Path(c.path_mol2).resolve()))
                frcmod_paths.append(str(Path(c.path_frcmod).resolve()))

        lig_resname = s.pose.res_name if s.pose else "LIG"
        plan_list.append(
            {
                "id": s.name,
                "rec_pdb": str(rec_pdb),
                "lig_sdf": lig_sdf,
                "lig_resname": lig_resname,
                "extra_sdfs": extra_sdfs,
                "extra_resnames": extra_resnames,
                "mol2_paths": mol2_paths,
                "frcmod_paths": frcmod_paths,
            }
        )

    plan_file = cwd / "md_setup_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_list, f, indent=2)

    console.print(
        f"[bold cyan]Dispatching topology setup for {len(plan_list)} System(s)...[/bold cyan]"
    )

    nf_args = build_nextflow_args(configs, cwd, entrypoint="md-setup")
    nf_args.extend(["--md_plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"
    ret_code = execute_nextflow(
        "md-setup", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print("[bold red]Nextflow MD setup failed.[/bold red]")
        sys.exit(ret_code)

    receipt_file = cwd / "md_setup_receipt.tsv"
    if not receipt_file.exists():
        console.print("[bold red]FATAL: md_setup_receipt.tsv not found.[/bold red]")
        sys.exit(1)

    sys_map = {s.name: s for s in raw_systems}
    registered = 0
    lines = receipt_file.read_text().strip().splitlines()
    for line in lines:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) == 5:
            sys_id, gro, top, ndx, posre_str = parts
            posre_files = posre_str.split(",") if posre_str else []
        elif len(parts) == 4:
            sys_id, gro, top, ndx = parts
            posre_files = []
        else:
            continue

        if sys_id in sys_map:
            s = sys_map[sys_id]
            s.path_gro = gro
            s.path_top = top
            s.path_ndx = ndx
            
            env = s.environment_dict
            env["posre_files"] = posre_files
            s.environment_dict = env
            
            s.state = "ready"
            registered += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ MD Setup completed! {registered} System(s) transitioned to 'ready'.[/bold green]"
    )
    console.print(
        "[dim]Topologies compiled. Run 'sbdd md' to launch production dynamics.[/dim]"
    )
