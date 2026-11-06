"""
File: cli/commands/rescore.py
Description: CLI handler for 'sbdd rescore': standalone Prime MM-GBSA and Consensus rescoring.
"""

import csv
import json
import sys
from pathlib import Path
from typing import Final
from rich.console import Console
from rich.table import Table
from cli.core.db import (
    DockedPose,
    DockingRun,
    ReceptorFrame,
    Substance,
    System,
    init_db,
)
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console: Final[Console] = Console()


def run_rescore(
    use_whitelist: bool,
    urgent: bool,
    local: bool,
    resume: bool,
) -> None:
    """Runs standalone Prime MM-GBSA and Consensus rescoring on the latest 3NF DockingRun."""
    cwd = Path.cwd()
    session = init_db()

    dock_run = (
        session.query(DockingRun)
        .order_by(DockingRun.created_at.desc())
        .first()
    )
    if not dock_run:
        session.close()
        console.print(
            "[bold red]FATAL: No DockingRun entities found in 3NF database. Run 'sbdd dock' or 'sbdd redock' first.[/bold red]"
        )
        sys.exit(1)

    frame_ids = sorted(
        list({dp.receptor_frame_id for dp in dock_run.docked_poses})
    )
    if not frame_ids:
        session.close()
        console.print(
            f"[bold red]FATAL: DockingRun '{dock_run.name}' contains zero DockedPose records.[/bold red]"
        )
        sys.exit(1)

    tasks: list[dict[str, str]] = []
    dock_summary_files: set[str] = set()

    for dp in dock_run.docked_poses:
        p_sdf = Path(dp.path_sdf).resolve()
        if not p_sdf.exists():
            session.close()
            raise FileNotFoundError(
                f"FATAL: Registered DockedPose SDF missing on disk: {p_sdf}"
            )
        dock_summary_files.add(str(p_sdf))

    for f_id in frame_ids:
        frame = session.query(ReceptorFrame).get(f_id)
        if not frame:
            continue
        frame_poses = [
            dp for dp in dock_run.docked_poses if dp.receptor_frame_id == frame.id
        ]
        if not frame_poses:
            continue

        frame_dir = Path(frame_poses[0].path_sdf).resolve().parent
        pv_file = frame_dir / f"{frame.name}_dock_pv.maegz"
        if not pv_file.exists() or pv_file.stat().st_size == 0:
            session.close()
            raise FileNotFoundError(
                f"FATAL: Glide Pose Viewer file '{pv_file}' not found for frame '{frame.name}' in DockingRun '{dock_run.name}'. "
                "Prime MM-GBSA rescoring requires Glide Pose Viewer (_dock_pv.maegz) output."
            )

        tasks.append({"name": frame.name, "pv_path": str(pv_file), "frame": frame.name})

        glide_csv = frame_dir / f"{frame.name}_dock.csv"
        if glide_csv.exists():
            dock_summary_files.add(str(glide_csv))

        if dock_run.mode == "self":
            for eng in ("glide", "unidock"):
                rmsd_csv = frame_dir / f"{frame.name}_{eng}_rmsd.csv"
                if rmsd_csv.exists():
                    dock_summary_files.add(str(rmsd_csv))

    if not tasks:
        session.close()
        console.print(
            "[bold red]FATAL: No Pose Viewer tasks constructed for Prime MM-GBSA rescoring.[/bold red]"
        )
        sys.exit(1)

    native_ligand_name = ""
    for s in session.query(System).all():
        if s.pose and s.pose.origin == "crystal" and s.pose.substance:
            native_ligand_name = s.pose.substance.name
            break

    dock_plan_file = cwd / "dock_plan.json"
    if dock_plan_file.exists():
        dock_summary_files.add(str(dock_plan_file.resolve()))

    whitelist_csv = cwd / "funnel_whitelist.csv"
    plan = {
        "use_whitelist": use_whitelist,
        "whitelist_path": (
            str(whitelist_csv.resolve()) if whitelist_csv.exists() else ""
        ),
        "dock_summary_files": sorted(list(dock_summary_files)),
        "tasks": tasks,
    }

    plan_file = cwd / "rescore_plan.json"
    with open(plan_file, "w", encoding="utf-8") as f:
        json.dump(plan, f, indent=2)

    console.print(
        f"[bold cyan]Generated rescore_plan.json for DockingRun '{dock_run.name}' ({len(tasks)} frame(s))...[/bold cyan]"
    )

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="rescore")
    nf_args.extend(["--rescore_plan", str(plan_file.resolve())])

    if native_ligand_name:
        if "--native_ligand" in nf_args:
            idx = nf_args.index("--native_ligand")
            nf_args[idx + 1] = native_ligand_name
        else:
            nf_args.extend(["--native_ligand", native_ligand_name])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    profile_name = "local" if local else "slurm"

    ret_code = execute_nextflow(
        "rescore", nf_args, cwd, profile=profile_name, resume=resume, configs=configs
    )

    if ret_code != 0:
        session.close()
        console.print("[bold red]Nextflow MM-GBSA rescoring failed.[/bold red]")
        sys.exit(ret_code)

    receipt_file = cwd / "rescore_receipt.tsv"
    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not generated by orchestrator. "
            "Database update aborted."
        )

    res_csv = cwd / "reports" / "ranking.csv"
    if not res_csv.exists() or res_csv.stat().st_size == 0:
        session.close()
        raise FileNotFoundError(
            f"FATAL: Consolidated ranking {res_csv} is missing. Process aborted."
        )

    with open(res_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw_title = row.get("Title")
            raw_frame = (
                row.get("Receptor")
                if row.get("Receptor") is not None
                else row.get("Frame")
            )
            title = raw_title.strip() if raw_title else ""
            frame_name = raw_frame.strip() if raw_frame else ""
            if not title or not frame_name:
                continue

            mmgbsa_raw = row.get("MMGBSA")
            strain_raw = row.get("Lig_Strain")
            ult_raw = row.get("Ultimate_Score")

            mmgbsa_val = (
                float(mmgbsa_raw)
                if (mmgbsa_raw and mmgbsa_raw.lower() != "nan")
                else None
            )
            strain_val = (
                float(strain_raw)
                if (strain_raw and strain_raw.lower() != "nan")
                else None
            )
            ult_val = (
                float(ult_raw)
                if (ult_raw and ult_raw.lower() != "nan")
                else None
            )

            r_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
            sub_obj = session.query(Substance).filter_by(name=title).first()
            if r_frame and sub_obj:
                matched_dps = (
                    session.query(DockedPose)
                    .filter_by(
                        docking_run_id=dock_run.id,
                        receptor_frame_id=r_frame.id,
                        substance_id=sub_obj.id,
                    )
                    .all()
                )
                for dp in matched_dps:
                    if mmgbsa_val is not None:
                        dp.mmgbsa_dg_bind = mmgbsa_val
                    if strain_val is not None:
                        dp.lig_strain_energy = strain_val
                    if ult_val is not None:
                        dp.consensus_score = ult_val

    session.commit()
    session.close()

    # Render summary results table
    table = Table(title="🎯 Post-Rescore Consensus Ranking Summary", show_header=True)
    table.add_column("Receptor", style="cyan")
    table.add_column("Ligand", style="bold")
    table.add_column("GlideScore", justify="right")
    table.add_column("UniDock Vina", justify="right")
    table.add_column("CNN Affinity", justify="right")
    table.add_column("MM-GBSA (kcal/mol)", justify="right", style="bold magenta")
    table.add_column("Glide RMSD (Å)", justify="center")
    table.add_column("UniDock RMSD (Å)", justify="center")
    table.add_column("Verdict", justify="center")

    with open(res_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rec = row.get("Receptor", row.get("Frame", "-"))
            lig = row.get("Title", "-")
            gscore = row.get("GlideScore", "-")
            uscore = row.get("UniDock_Vina", "-")
            cnn = row.get("Glide_CNNaffinity") or row.get("UniDock_CNNaffinity") or "-"
            mmgbsa = row.get("MMGBSA", "-")
            g_rmsd = row.get("InSitu_RMSD_Glide_A", "-")
            n_rmsd = row.get("InSitu_RMSD_UniDock_A", "-")

            def _fmt(val):
                try:
                    return f"{float(val):.2f}"
                except (ValueError, TypeError):
                    return str(val) if val and val != "nan" else "-"

            def _color_rmsd(val_str):
                if val_str and val_str not in ["nan", "-"]:
                    try:
                        v = float(val_str)
                        return (
                            f"[bold green]{v:.2f}[/bold green]"
                            if v <= 2.0
                            else f"[bold red]{v:.2f}[/bold red]"
                        )
                    except ValueError:
                        return val_str
                return "-"

            verdict = "[dim]-[/dim]"
            try:
                vg = float(g_rmsd) if g_rmsd and g_rmsd not in ["nan", "-"] else None
                vn = float(n_rmsd) if n_rmsd and n_rmsd not in ["nan", "-"] else None
                if (vg is not None and vg <= 2.5) or (vn is not None and vn <= 2.5):
                    verdict = "[bold green]✓ SUCCESS[/bold green]"
                elif vg is not None or vn is not None:
                    verdict = "[bold yellow]⚠ DRIFT[/bold yellow]"
            except ValueError:
                pass

            table.add_row(
                rec,
                lig,
                _fmt(gscore),
                _fmt(uscore),
                _fmt(cnn),
                _fmt(mmgbsa),
                _color_rmsd(g_rmsd),
                _color_rmsd(n_rmsd),
                verdict,
            )

    console.print(table)
    console.print(
        "\n[bold green]✓ Prime MM-GBSA rescoring completed and synchronized with SQLite database![/bold green]\n"
    )
