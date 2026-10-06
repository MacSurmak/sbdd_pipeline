"""
File: cli/commands/dock.py
Description: CLI handler for 'sbdd dock'. Reads DB to construct Cross-Docking plan.
"""

import json
import sys
from pathlib import Path
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_dock(args):
    cwd = Path.cwd()
    session = init_db()
    is_self = getattr(args, "self_dock", False)

    # 1. Fetch requested Ready Frames
    query = session.query(ReceptorFrame).filter_by(state="ready")
    if getattr(args, "target", None):
        target_obj = session.query(Target).filter_by(name=args.target).first()
        if not target_obj:
            console.print(
                f"[bold red]FATAL: Target '{args.target}' not found.[/bold red]"
            )
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif getattr(args, "frames", None) or getattr(args, "receptors", None):
        req = getattr(args, "frames", None) or getattr(args, "receptors", None)
        f_names = [f.strip() for f in req.split(",") if f.strip()]
        query = query.filter(ReceptorFrame.name.in_(f_names))

    frames = query.all()
    if not frames:
        console.print(
            "[bold red]FATAL: No ready ReceptorFrames found to dock. Run 'sbdd prep' first![/bold red]"
        )
        sys.exit(1)

    # 2. Build Plan (Self 1-to-1 vs Matrix Screening)
    plan = {
        "mode": "self" if is_self else "matrix",
        "receptors": [],
        "ligands": [],
        "pairs": [],
    }

    if is_self:
        active_frames = []

        for r in frames:
            target_sys = (
                session.query(System)
                .filter(System.receptor_frame_id == r.id, System.pose_id.isnot(None))
                .first()
            )
            if not target_sys or not target_sys.pose:
                continue

            native_pose = target_sys.pose
            native_lig = native_pose.ligand

            active_frames.append((r, native_pose, native_lig))

        if not active_frames:
            console.print(
                f"[bold green]✓ All {len(frames)} ready frame(s) already have redocking poses completed. Nothing to do![/bold green]"
            )
            session.close()
            return

        for r, native_pose, native_lig in active_frames:
            # Reconcile any existing corrupted scratch paths strictly with the permanent prepared/ folder
            frame_p = Path(r.path).resolve()
            if "/work/" in str(frame_p) or not frame_p.exists():
                canonical_mae = (cwd / "prepared" / frame_p.name).resolve()
                if canonical_mae.exists():
                    frame_p = canonical_mae
                    r.path = str(canonical_mae)
                    session.add(r)
                else:
                    raise FileNotFoundError(
                        f"FATAL: Receptor frame '{r.name}' references invalid path: {r.path}. "
                        f"File not found at canonical location: {canonical_mae}"
                    )

            path_mae = frame_p.with_suffix(".mae")
            path_pdb = frame_p.with_suffix(".pdb")

            if not path_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Receptor MAE not found for frame '{r.name}': {path_mae}"
                )
            if not path_pdb.exists():
                raise FileNotFoundError(
                    f"FATAL: Receptor PDB not found for frame '{r.name}': {path_pdb}"
                )

            # Reconcile pose scratch paths strictly with the permanent prepared/ folder
            pose_p = Path(native_pose.path).resolve()
            if "/work/" in str(pose_p) or not pose_p.exists():
                canonical_sdf = (cwd / "prepared" / pose_p.name).resolve()
                if canonical_sdf.exists():
                    pose_p = canonical_sdf
                    native_pose.path = str(canonical_sdf)
                    session.add(native_pose)
                else:
                    raise FileNotFoundError(
                        f"FATAL: Pose '{native_pose.name}' references invalid path: {native_pose.path}. "
                        f"File not found at canonical location: {canonical_sdf}"
                    )

            pose_sdf = pose_p.with_suffix(".sdf")
            pose_mae = pose_p.with_suffix(".mae")

            if not pose_sdf.exists():
                raise FileNotFoundError(
                    f"FATAL: Reference pose SDF not found for frame '{r.name}': {pose_sdf}"
                )
            if not pose_mae.exists():
                raise FileNotFoundError(
                    f"FATAL: Reference pose MAE not found for frame '{r.name}': {pose_mae}"
                )

            plan["receptors"].append(
                {
                    "name": r.name,
                    "path_mae": str(path_mae),
                    "path_pdb": str(path_pdb),
                    "ref_sdf": str(pose_sdf),
                    "ref_mae": str(pose_mae),
                    "native_ligand": native_lig.name,
                }
            )

            unique_lig_name = f"{r.name}_{native_lig.name}"

            # Strict validation: Candidate ligand MUST be in state 'ready' (processed via LigPrep)
            if (
                native_lig.state != "ready"
                or not native_lig.path_mae
                or not native_lig.path_sdf
            ):
                console.print(
                    f"[bold red]FATAL: Chemical Ligand '{native_lig.name}' is in state '{native_lig.state}'! "
                    f"Run 'sbdd prep' first so LigPrep generates canonical 3D states.[/bold red]"
                )
                session.close()
                sys.exit(1)

            lig_mae = Path(native_lig.path_mae).resolve()
            lig_sdf = Path(native_lig.path_sdf).resolve()

            if not lig_mae.exists() or not lig_sdf.exists():
                console.print(
                    f"[bold red]FATAL: LigPrep deliverables missing for '{native_lig.name}': {lig_mae}. Run 'sbdd prep'.[/bold red]"
                )
                session.close()
                sys.exit(1)

            plan["ligands"].append(
                {
                    "name": unique_lig_name,
                    "path_sdf": str(lig_sdf),
                    "path_mae": str(lig_mae),
                    "is_native": True,
                    "receptor": r.name,
                    "ligand_title": native_lig.name,
                }
            )

        console.print(
            f"[bold green]✓ 1-to-1 Redocking Mode: Matched {len(active_frames)} receptors directly to their native ligands.[/bold green]"
        )

    else:
        if not getattr(args, "ligands", None):
            console.print(
                "[bold red]FATAL: Parameter --ligands is required for cross-docking matrix mode.[/bold red]"
            )
            sys.exit(1)

        lig_names = [n.strip() for n in args.ligands.split(",")]
        ligands = session.query(Ligand).filter(Ligand.name.in_(lig_names)).all()
        if len(ligands) != len(lig_names):
            found = {l.name for l in ligands}
            missing = set(lig_names) - found
            console.print(
                f"[bold red]FATAL: Ligands not found in DB: {missing}[/bold red]"
            )
            sys.exit(1)

        for r in frames:
            sys_obj = session.query(System).filter_by(receptor_frame_id=r.id).first()
            if not sys_obj or not sys_obj.pose:
                console.print(
                    f"[bold red]FATAL: Receptor '{r.name}' has no reference pose! Screening requires a pocket reference.[/bold red]"
                )
                sys.exit(1)

            pose_obj = sys_obj.pose
            native_lig = pose_obj.ligand.name
            pose_p = Path(pose_obj.path).resolve()

            ref_sdf = pose_p.with_suffix(".sdf")
            ref_mae = pose_p.with_suffix(".mae")

            if not ref_sdf.exists() or not ref_mae.exists():
                console.print(
                    f"[bold red]FATAL: Reference pose files missing for '{r.name}'. Need both .sdf and .mae![/bold red]"
                )
                sys.exit(1)

            frame_p = Path(r.path).resolve()
            path_mae = frame_p.with_suffix(".mae")
            path_pdb = frame_p.with_suffix(".pdb")

            if not path_mae.exists() or not path_pdb.exists():
                console.print(
                    f"[bold red]FATAL: Receptor files missing for '{r.name}'. Need both .pdb and .mae![/bold red]"
                )
                sys.exit(1)

            plan["receptors"].append(
                {
                    "name": r.name,
                    "path_mae": str(path_mae),
                    "path_pdb": str(path_pdb),
                    "ref_sdf": str(ref_sdf),
                    "ref_mae": str(ref_mae),
                    "native_ligand": native_lig,
                }
            )

        for l in ligands:
            if l.state != "ready" or not l.path_mae or not l.path_sdf:
                console.print(
                    f"[bold red]FATAL: Chemical Ligand '{l.name}' is in state '{l.state}'! Run 'sbdd prep' first.[/bold red]"
                )
                session.close()
                sys.exit(1)

            lig_mae = Path(l.path_mae).resolve()
            lig_sdf = Path(l.path_sdf).resolve()

            if not lig_mae.exists() or not lig_sdf.exists():
                console.print(
                    f"[bold red]FATAL: Prepared files missing for '{l.name}' ({lig_mae}). Run 'sbdd prep' first.[/bold red]"
                )
                session.close()
                sys.exit(1)

            has_crystal_pose = (
                session.query(Pose).filter_by(ligand_id=l.id, origin="crystal").first()
                is not None
            )

            plan["ligands"].append(
                {
                    "name": l.name,
                    "path_sdf": str(lig_sdf),
                    "path_mae": str(lig_mae),
                    "is_native": has_crystal_pose,
                }
            )

        console.print(
            f"[bold cyan]Cross-Docking Matrix Plan generated: {len(frames)} Receptors × {len(ligands)} Ligand Sets[/bold cyan]"
        )

    plan_file = cwd / "dock_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan, f, indent=2)

    # 4. Build Nextflow args
    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="screen")
    nf_args.extend(["--dock_plan", str(plan_file.resolve())])

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

    # 6. Register results in Database strictly via Plan & Receipt protocol
    num_ligs = len(plan.get("ligands", []))
    default_name = (
        f"redock_{len(frames)}" if is_self else f"screen_{len(frames)}x{num_ligs}"
    )
    poses_name = args.name if getattr(args, "name", None) else default_name

    receipt_file = cwd / "dock_receipt.tsv"
    if not receipt_file.exists() or receipt_file.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Execution receipt {receipt_file} was not generated by Nextflow. "
            "Docking failed or terminated without producing deliverables."
        )

    res_csv = cwd / "reports" / "ranking.csv"
    if not res_csv.exists() or res_csv.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Ranking summary CSV {res_csv} does not exist. Pipeline failed to summarize docking."
        )

    # Parse receipt deterministically: map of (frame_id, engine) -> absolute SDF path
    receipt_map = {}
    with open(receipt_file, "r", encoding="utf-8") as rf:
        for line in rf:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                r_frame_id, r_engine, r_sdf_path = parts[0], parts[1], parts[2]
                receipt_map[(r_frame_id, r_engine)] = r_sdf_path

    import csv

    poses_count = 0
    with open(res_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("Title", "").strip()
            frame_name = row.get("Receptor", row.get("Frame", "")).strip()
            if not title or not frame_name:
                raise ValueError(
                    f"FATAL: Corrupted ranking row with missing Title or Frame: {row}"
                )

            r_frame = session.query(ReceptorFrame).filter_by(name=frame_name).first()
            if not r_frame:
                raise ValueError(
                    f"FATAL: ReceptorFrame '{frame_name}' in ranking is not registered in SQLite database!"
                )

            lig = session.query(Ligand).filter_by(name=title).first()
            if not lig:
                raise ValueError(
                    f"FATAL: Ligand '{title}' in ranking is not registered in SQLite database!"
                )

            # Resolve canonical SDF path strictly via the orchestrator receipt
            canonical_sdf = (
                receipt_map.get((frame_name, "glide"))
                or receipt_map.get((frame_name, "gnina"))
            )
            if not canonical_sdf:
                raise FileNotFoundError(
                    f"FATAL: No receipt record found in dock_receipt.tsv for receptor frame '{frame_name}'!"
                )

            sdf_path = Path(canonical_sdf).resolve()
            if not sdf_path.exists() or sdf_path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Pose deliverable specified by receipt does not exist: {sdf_path}"
                )

            pose_name = f"{frame_name}_{title}_dock"
            scores_dict = {
                k: float(v)
                for k, v in row.items()
                if k not in ["Rank", "Title", "Receptor", "Frame"]
                and v
                and v != "nan"
                and v != ""
            }

            ex = session.query(Pose).filter_by(name=pose_name).first()
            if ex:
                ex.path = str(sdf_path)
                ex.scores_dict = scores_dict
                ex.state = "ready"
            else:
                new_p = Pose(
                    ligand_id=lig.id,
                    receptor_frame_id=r_frame.id,
                    name=pose_name,
                    origin="docking",
                    state="ready",
                    path=str(sdf_path),
                )
                new_p.scores_dict = scores_dict
                session.add(new_p)
            poses_count += 1

    session.commit()
    console.print(
        f"\n[bold green]✓ Docking verified via receipt! Registered {poses_count} pose(s) in DB as '{poses_name}'.[/bold green]\n"
    )

        # Render Rich Summary Table
        try:
            from rich.table import Table

            table = Table(title="🎯 Redocking & Scoring Summary", show_header=True)
            table.add_column("Receptor", style="cyan")
            table.add_column("Ligand", style="bold")
            table.add_column("GlideScore", justify="right")
            table.add_column("CNN Affinity", justify="right")
            table.add_column("MM-GBSA (kcal/mol)", justify="right")
            table.add_column("Glide RMSD (Å)", justify="center")
            table.add_column("GNINA RMSD (Å)", justify="center")
            table.add_column("Verdict", justify="center")

            with open(res_csv, "r", encoding="utf-8", errors="ignore") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rec = row.get("Receptor", row.get("Frame", "-"))
                    lig = row.get("Title", "-")
                    gscore = row.get("GlideScore", "-")
                    cnn = (
                        row.get("Resc_CNNaffinity")
                        or row.get("Dock_CNNaffinity")
                        or "-"
                    )
                    mmgbsa = row.get("MMGBSA", "-")
                    g_rmsd = row.get("InSitu_RMSD_Glide_A", "-")
                    n_rmsd = row.get("InSitu_RMSD_GNINA_A", "-")

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

                    disp_glide = _color_rmsd(g_rmsd)
                    disp_gnina = _color_rmsd(n_rmsd)

                    # Determine verdict: SUCCESS if either engine reached <= 2.5 A (allosteric tolerance)
                    verdict = "[dim]-[/dim]"
                    try:
                        vg = (
                            float(g_rmsd)
                            if g_rmsd and g_rmsd not in ["nan", "-"]
                            else None
                        )
                        vn = (
                            float(n_rmsd)
                            if n_rmsd and n_rmsd not in ["nan", "-"]
                            else None
                        )
                        if (vg is not None and vg <= 2.5) or (
                            vn is not None and vn <= 2.5
                        ):
                            verdict = "[bold green]✓ SUCCESS[/bold green]"
                        elif vg is not None or vn is not None:
                            verdict = "[bold yellow]⚠ DRIFT[/bold yellow]"
                    except ValueError:
                        pass

                    table.add_row(
                        rec,
                        lig,
                        _fmt(gscore),
                        _fmt(cnn),
                        _fmt(mmgbsa),
                        disp_glide,
                        disp_gnina,
                        verdict,
                    )
            console.print(table)
        except Exception as e:
            console.print(
                f"[yellow]Warning: Could not render console table: {e}[/yellow]"
            )

    session.close()
