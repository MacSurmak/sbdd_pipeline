"""
File: cli/commands/prep.py
Description: CLI handler for 'sbdd prep'. Reads DB, builds Execution Plan, runs NF, and updates DB.
"""

import json
from pathlib import Path
import sys
from rich.console import Console
from cli.core.db import init_db, Receptor
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_prep(args):
    cwd = Path.cwd()
    session = init_db()

    # 1. Fetch receptors needing preparation
    query = session.query(Receptor).filter(Receptor.status == "unprepared")
    if getattr(args, "receptors", None):
        target_names = [n.strip() for n in args.receptors.split(",")]
        query = query.filter(Receptor.name.in_(target_names))

    receptors = query.all()
    if not receptors:
        console.print(
            "[green]No unprepared receptors found. Everything is up to date![/green]"
        )
        session.close()
        return

    # 2. Build Execution Plan
    plan_data = []
    for r in receptors:
        plan_data.append(
            {
                "name": r.name,
                "path": r.path,
                "prep_mode": r.prep_mode,
                "cofactors": r.meta_dict.get("cofactors", []),
                "native_ligand": r.meta_dict.get("native_ligand", ""),
            }
        )

    plan_file = cwd / "prep_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_data, f, indent=2)

    console.print(
        f"[bold cyan]Generated prep_plan.json for {len(receptors)} receptors.[/bold cyan]"
    )

    # 3. Build Nextflow args and execute
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="prepare")
    nf_args.extend(["--plan", str(plan_file.resolve())])

    if getattr(args, "urgent", False):
        nf_args.extend(["--queue", "urgent"])

    profile = "local" if getattr(args, "local", False) else "slurm"
    resume = not getattr(args, "no_resume", False)

    ret_code = execute_nextflow(
        "prepare", nf_args, cwd, profile=profile, resume=resume, configs=configs
    )

    if ret_code != 0:
        console.print(
            "[bold red]Nextflow execution failed. Database not updated.[/bold red]"
        )
        sys.exit(ret_code)

    # 4. Read Receipt and update Database
    receipt_file = cwd / "prep_receipt.tsv"
    if not receipt_file.exists():
        console.print(
            "[bold red]FATAL: prep_receipt.tsv not found! Nextflow might have failed silently.[/bold red]"
        )
        sys.exit(1)

    lines = receipt_file.read_text().strip().splitlines()
    success_count = 0

    # DB models needed for auto-registering native ligands
    from cli.core.db import Ligand

    for line in lines:
        if not line.strip():
            continue
        name, raw_out_path = line.split("\t")

        # Resolve the published path (in prepared/ folder) instead of work/ cache
        published_path = cwd / "prepared" / Path(raw_out_path).name
        if not published_path.exists():
            published_path = Path(raw_out_path)  # Fallback

        ready_name = f"{name}_ready"
        existing = session.query(Receptor).filter_by(name=ready_name).first()

        if existing:
            existing.path = str(published_path)
            existing.status = "ready"
        else:
            new_rec = Receptor(
                name=ready_name,
                path=str(published_path),
                status="ready",
                prep_mode="ready",
                parent_name=name,
            )
            session.add(new_rec)

        orig = session.query(Receptor).filter_by(name=name).first()
        if orig:
            orig.status = "processed"

            # Auto-register extracted native ligand if it exists
            native_lig = orig.meta_dict.get("native_ligand")
            if native_lig:
                lig_path = cwd / "prepared" / f"{name}_ref_native_prepared.sdf"
                if not lig_path.exists():
                    lig_path = cwd / "prepared" / f"{name}_ref_native.sdf"

                if lig_path.exists():
                    lig_name = f"{name}_native"
                    ex_lig = session.query(Ligand).filter_by(name=lig_name).first()
                    if not ex_lig:
                        new_lig = Ligand(
                            name=lig_name,
                            path=str(lig_path),
                            status="ready",
                            is_native=True,
                        )
                        session.add(new_lig)

        success_count += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Successfully prepared {success_count} receptors and updated the database![/bold green]"
    )
    console.print("[dim]Use 'sbdd tree' to view the updated lineage.[/dim]")
