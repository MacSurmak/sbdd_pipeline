"""
File: cli/commands/prep.py
Description: CLI handler for 'sbdd prep' using clean 3NF database lineage.
Prepares raw ReceptorFrames and in-situ Poses into ready states without mutating raw records.
"""

import json
from pathlib import Path
import sys
from rich.console import Console
from cli.core.db import init_db, Target, ReceptorFrame, Ligand, Pose, System
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def cmd_prep(args):
    cwd = Path.cwd()
    session = init_db()

    # 1. Fetch raw frames that need preparation
    query = session.query(ReceptorFrame).filter(ReceptorFrame.state == "raw")
    if getattr(args, "target", None):
        target_obj = session.query(Target).filter_by(name=args.target).first()
        if not target_obj:
            console.print(
                f"[bold red]FATAL: Target '{args.target}' not found in database.[/bold red]"
            )
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif getattr(args, "frames", None):
        frame_names = [f.strip() for f in args.frames.split(",") if f.strip()]
        query = query.filter(ReceptorFrame.name.in_(frame_names))

    raw_frames = query.all()

    # 2. Fetch raw chemical ligands that need LigPrep
    query_ligands = session.query(Ligand).filter_by(state="raw").all()

    if not raw_frames and not query_ligands:
        console.print(
            "[green]No raw frames or ligands needing preparation found. Everything is up to date![/green]"
        )
        session.close()
        return

    prepared_dir = (cwd / "prepared").resolve()
    prepared_lig_dir = prepared_dir / "ligands"
    raw_lig_dir = prepared_lig_dir / "raw"
    raw_lig_dir.mkdir(parents=True, exist_ok=True)

    lig_plan = []
    for lig in query_ligands:
        lig_raw_sdf = raw_lig_dir / f"{lig.name}.sdf"
        # Always export 2D/3D seed SDF from canonical SMILES for un-biased LigPrep
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            m = Chem.MolFromSmiles(lig.smiles)
            if m:
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m, randomSeed=42)
                m.SetProp("_Name", lig.name)
                w = Chem.SDWriter(str(lig_raw_sdf))
                w.write(m)
                w.close()
        except Exception:
            pass

        if lig_raw_sdf.exists() and lig_raw_sdf.stat().st_size > 0:
            lig_plan.append({
                "name": lig.name,
                "path": str(lig_raw_sdf.resolve()),
                "status": "unprepared"
            })

    if not raw_frames and not lig_plan:
        console.print(
            "[green]No raw frames or ligands needing preparation found. Everything is up to date![/green]"
        )
        session.close()
        return

    # 3. Build Plan for Nextflow
    rec_plan = []
    frame_map = {}

    for frame in raw_frames:
        frame_map[frame.name] = frame
        system = session.query(System).filter_by(receptor_frame_id=frame.id).first()
        native_lig_code = ""
        if system and system.pose:
            native_lig_code = system.pose.ligand.name

        prep_mode = "light" if frame.source_type == "md_centroid" else "full"
        rec_plan.append(
            {
                "name": frame.name,
                "path": str(Path(frame.path).resolve()),
                "prep_mode": prep_mode,
                "cofactors": frame.cofactors_list,
                "native_ligand": native_lig_code,
            }
        )

    plan_data = {
        "receptors": rec_plan,
        "ligands": lig_plan
    }

    plan_file = cwd / "prep_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_data, f, indent=2)

    console.print(
        f"[bold cyan]Generated prep_plan.json ({len(rec_plan)} receptor frame(s), {len(lig_plan)} ligand(s))...[/bold cyan]"
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
            "[bold red]Nextflow preparation failed. Database not updated.[/bold red]"
        )
        sys.exit(ret_code)

    # 4. Read Receipts and update database records
    receipt_file = cwd / "prep_receipt.tsv"
    prepared_count = 0

    if receipt_file.exists():
        lines = receipt_file.read_text().strip().splitlines()
        for line in lines:
            if not line.strip():
                continue
            orig_frame_name, mae_filename = line.split("\t")
            frame = frame_map.get(orig_frame_name)
            if not frame:
                continue

            prep_mae_path = (cwd / "prepared" / Path(mae_filename).name).resolve()
            if not prep_mae_path.exists():
                raise FileNotFoundError(
                    f"FATAL: Prepared receptor MAE not found at canonical location: {prep_mae_path}"
                )

            frame.path = str(prep_mae_path)
            frame.state = "ready"

            # In-situ reference pose from PrepWizard
            system = session.query(System).filter_by(receptor_frame_id=frame.id).first()
            if system and system.pose:
                prep_sdf_path = (
                    cwd / "prepared" / f"{orig_frame_name}_ref_native_prepared.sdf"
                ).resolve()
                if not prep_sdf_path.exists():
                    prep_sdf_path = (
                        cwd / "prepared" / f"{orig_frame_name}_ref_native.sdf"
                    ).resolve()

                if prep_sdf_path.exists() and prep_sdf_path.stat().st_size > 0:
                    system.pose.path = str(prep_sdf_path)
                    system.pose.state = "ready"

                    # Calculate actual formal charge from the prepared in-situ ligand
                    try:
                        from rdkit import Chem
                        suppl = Chem.SDMolSupplier(str(prep_sdf_path), removeHs=False)
                        mol = next((m for m in suppl if m is not None), None)
                        if mol:
                            ch = Chem.GetFormalCharge(mol)
                            system.pose.scores_dict = {"formal_charge": ch}
                    except Exception:
                        pass

            if system:
                system.state = "ready"
            prepared_count += 1

    # 5. Update prepared Ligands in DB with LigPrep deliverable paths
    lig_receipt = cwd / "prep_ligand_receipt.tsv"
    prep_lig_count = 0
    if lig_receipt.exists():
        for line in lig_receipt.read_text().strip().splitlines():
            if not line.strip():
                continue
            lig_name, mae_name = line.split("\t")
            lig_obj = session.query(Ligand).filter_by(name=lig_name).first()
            if lig_obj:
                prep_mae = (cwd / "ligands" / mae_name).resolve()
                sdf_name = mae_name.replace(".maegz", ".sdf").replace(".mae", ".sdf")
                prep_sdf = (cwd / "ligands" / sdf_name).resolve()

                lig_obj.state = "ready"
                lig_obj.path_mae = str(prep_mae)
                lig_obj.path_sdf = str(prep_sdf)
                prep_lig_count += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Successfully prepared {prepared_count} receptor frame(s) and {prep_lig_count} chemical ligand(s)![/bold green]"
    )
    console.print(
        "[dim]Crystal in-situ poses and LigPrep 3D states are ready. Use 'sbdd tree' to view lineage.[/dim]"
    )
