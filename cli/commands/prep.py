"""
File: cli/commands/prep.py
Description: CLI handler for 'sbdd prep' using clean 3NF database lineage.
Prepares raw ReceptorFrames and in-situ Poses into ready states without mutating raw records.
"""

import json
from pathlib import Path
import sys
from rich.console import Console
from cli.core.db import (
    init_db,
    Target,
    ReceptorFrame,
    Substance,
    SubstanceState,
    Pose,
    System,
)
from cli.core.runner import execute_nextflow
from cli.core.config_loader import load_project_configs, build_nextflow_args

console = Console()


def run_prep(
    target: str | None = None,
    frames: list[str] | None = None,
    urgent: bool = False,
    local: bool = False,
    resume: bool = True,
) -> None:
    """Orchestrates Tri-partite preparation: Receptors (PrepWizard), Substances (LigPrep), and Poses (Antechamber)."""
    cwd = Path.cwd()
    session = init_db()

    # 1. Fetch Receptors needing PrepWizard (Strictly exclude already equilibrated MD centroids)
    query = session.query(ReceptorFrame).filter(
        ReceptorFrame.state == "raw",
        ReceptorFrame.source_type != "md_centroid",
    )
    if target:
        target_obj = session.query(Target).filter_by(name=target).first()
        if not target_obj:
            console.print(
                f"[bold red]FATAL: Target '{target}' not found in database.[/bold red]"
            )
            sys.exit(1)
        query = query.filter(ReceptorFrame.target_id == target_obj.id)
    elif frames:
        query = query.filter(ReceptorFrame.name.in_(frames))

    raw_frames = query.all()

    # 2. Fetch abstract library Substances needing LigPrep (only type = 'ligand')
    query_substances = (
        session.query(Substance).filter_by(state="raw", substance_type="ligand").all()
    )

    if not raw_frames and not query_substances:
        console.print(
            "[green]No raw frames or substances found to prepare. Everything is up to date![/green]"
        )
        session.close()
        return

    # Build Plan
    rec_plan = []
    frame_map = {}
    for frame in raw_frames:
        frame_map[frame.name] = frame
        system = session.query(System).filter_by(receptor_frame_id=frame.id).first()
        heterogens = []
        if system and system.pose:
            heterogens.append(
                {
                    "res_name": system.pose.res_name,
                    "pose_name": system.pose.name,
                    "inchikey": system.pose.substance.inchikey,
                    "is_ligand": True,
                    "is_ion": False,
                }
            )
        if system:
            for cof in system.cofactors:
                heterogens.append(
                    {
                        "res_name": cof.res_name,
                        "pose_name": cof.name,
                        "inchikey": cof.substance.inchikey,
                        "is_ligand": False,
                        "is_ion": (cof.substance.substance_type == "ion"),
                    }
                )

        native_lig_code = system.pose.res_name if (system and system.pose) else ""
        cofactor_names = [h["res_name"] for h in heterogens if not h["is_ligand"]]

        rec_plan.append(
            {
                "name": frame.name,
                "target_name": frame.target.name,
                "path": str(Path(frame.path).resolve()),
                "prep_mode": "full",
                "cofactors": cofactor_names,
                "native_ligand": native_lig_code,
                "heterogens": heterogens,
            }
        )

    sub_plan = []
    if query_substances:
        for sub in query_substances:
            if not sub.path_sdf:
                session.close()
                raise ValueError(
                    f"FATAL: Substance '{sub.name}' has no registered path_sdf. Database lineage is corrupted."
                )

            sub_sdf = Path(sub.path_sdf).resolve()
            if not sub_sdf.exists() or sub_sdf.stat().st_size == 0:
                session.close()
                raise FileNotFoundError(
                    f"FATAL: Registered reference SDF for substance '{sub.name}' is missing or empty on disk: {sub_sdf}"
                )

            sub_plan.append(
                {
                    "name": sub.name,
                    "inchikey": sub.inchikey,
                    "path": str(sub_sdf),
                    "status": "unprepared",
                }
            )

    plan_data = {"receptors": rec_plan, "ligands": sub_plan}
    plan_file = cwd / "prep_plan.json"
    with open(plan_file, "w") as f:
        json.dump(plan_data, f, indent=2)

    console.print(
        f"[bold cyan]Generated prep_plan.json ({len(rec_plan)} receptors, {len(sub_plan)} substances).[/bold cyan]"
    )

    # Execute Nextflow (Purge stale receipts to prevent cross-run pollution)
    (cwd / "prep_receipt.tsv").unlink(missing_ok=True)
    (cwd / "prep_ligand_receipt.tsv").unlink(missing_ok=True)

    configs = load_project_configs(cwd)
    nf_args = build_nextflow_args(configs, cwd, entrypoint="prepare")
    nf_args.extend(["--plan", str(plan_file.resolve())])

    if urgent:
        nf_args.extend(["--queue", "urgent"])

    ret_code = execute_nextflow(
        "prepare",
        nf_args,
        cwd,
        profile="local" if local else "slurm",
        resume=resume,
        configs=configs,
    )
    if ret_code != 0:
        console.print(
            "[bold red]Nextflow preparation failed. Database not updated.[/bold red]"
        )
        sys.exit(ret_code)

    # Read Receipts strictly for components actually included in this plan execution
    prep_rec_count, prep_sub_count = 0, 0
    rec_receipt = cwd / "prep_receipt.tsv"
    if rec_plan and rec_receipt.exists() and rec_receipt.stat().st_size > 0:
        for line in rec_receipt.read_text().strip().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(
                    f"FATAL: Corrupted line in {rec_receipt}: '{line}'. Expected 3 tab-separated fields."
                )
            orig_name, r_mae, r_pdb = parts
            frame = frame_map.get(orig_name)
            if not frame:
                raise KeyError(
                    f"FATAL: Prepared frame '{orig_name}' in receipt is not present in plan!"
                )

            mae_path = Path(r_mae).resolve()
            pdb_path = Path(r_pdb).resolve()
            if not mae_path.exists() or mae_path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Prepared receptor MAE missing on disk: {mae_path}"
                )
            if not pdb_path.exists() or pdb_path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Prepared receptor PDB missing on disk: {pdb_path}"
                )

            frame.path_mae = str(mae_path)
            frame.path_pdb = str(pdb_path)
            frame.path = str(pdb_path)
            frame.state = "ready"
            prep_rec_count += 1

    # Synchronize microstates and poses prepared by PrepWizard
    poses_receipt = cwd / "prep_poses_receipt.tsv"
    if rec_plan and poses_receipt.exists() and poses_receipt.stat().st_size > 0:
        from cli.core.chem import get_microstate_identity
        import shutil

        for line in poses_receipt.read_text().strip().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            p_name, inchikey, res_clean, charge_str, rel_sdf, rel_mae = parts
            sdf_path = (cwd / rel_sdf).resolve() if not Path(rel_sdf).is_absolute() else Path(rel_sdf).resolve()
            mae_path = (cwd / rel_mae).resolve() if not Path(rel_mae).is_absolute() else Path(rel_mae).resolve()

            if not sdf_path.exists() or sdf_path.stat().st_size == 0:
                continue

            (
                parent_inchikey,
                parent_smiles,
                mol_name,
                state_hash,
                isomeric_smiles,
                formal_charge,
            ) = get_microstate_identity(sdf_path)

            sub_obj = session.query(Substance).filter_by(inchikey=parent_inchikey).first()
            if not sub_obj:
                sub_obj = session.query(Substance).filter_by(inchikey=inchikey).first()

            if sub_obj:
                state_entity = (
                    session.query(SubstanceState)
                    .filter_by(substance_id=sub_obj.id, state_hash=state_hash)
                    .first()
                )
                if not state_entity:
                    state_dir = cwd / "substances" / sub_obj.inchikey / "states" / state_hash
                    state_dir.mkdir(parents=True, exist_ok=True)
                    state_sdf_file = state_dir / f"{res_clean}_state.sdf"
                    shutil.copy2(str(sdf_path), str(state_sdf_file))

                    state_entity = SubstanceState(
                        substance_id=sub_obj.id,
                        state_hash=state_hash,
                        isomeric_smiles=isomeric_smiles,
                        formal_charge=formal_charge,
                        state="ready",
                        path_sdf=str(state_sdf_file.resolve()),
                    )
                    session.add(state_entity)
                    session.flush()

                pose_obj = session.query(Pose).filter_by(name=p_name).first()
                if pose_obj:
                    pose_obj.substance_state_id = state_entity.id
                    pose_obj.path_sdf = str(sdf_path)
                    pose_obj.path = str(sdf_path)
                    if mae_path and mae_path.exists():
                        pose_obj.path_mae = str(mae_path)
                    pose_obj.state = "ready"

    sub_receipt = cwd / "prep_ligand_receipt.tsv"
    if sub_plan and sub_receipt.exists() and sub_receipt.stat().st_size > 0:
        for line in sub_receipt.read_text().strip().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(
                    f"FATAL: Corrupted line in {sub_receipt}: '{line}'. Expected 3 tab-separated fields."
                )
            sub_name, rel_mae, rel_sdf = parts
            sub_obj = session.query(Substance).filter_by(name=sub_name).first()
            if not sub_obj:
                raise KeyError(
                    f"FATAL: Prepared substance '{sub_name}' in receipt is not registered in DB!"
                )

            mae_path = Path(rel_mae).resolve()
            sdf_path = Path(rel_sdf).resolve()
            if not mae_path.exists() or mae_path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Prepared ligand MAE missing: {mae_path}"
                )
            if not sdf_path.exists() or sdf_path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Prepared ligand SDF missing: {sdf_path}"
                )

            sub_obj.path_mae = str(mae_path)
            sub_obj.path_sdf = str(sdf_path)
            sub_obj.state = "ready"
            prep_sub_count += 1

    session.commit()
    session.close()

    console.print(
        f"\n[bold green]✓ Verified deliverables via receipts: {prep_rec_count} receptor(s), "
        f"{prep_sub_count} substance(s) prepared![/bold green]"
    )
