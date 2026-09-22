"""
File: src/utils/lid_export.py
Description: Orchestrates 2D Ligand Interaction Diagram (LID) generation for
the Top-N screened hits, native reference ligand, and benchmark controls.
"""
import json
import os
from pathlib import Path
from typing import List, Dict, Any, Optional
import pandas as pd
from loguru import logger

from utils.cmd_runner import run_command


def _resolve_receptor_file(cwd: Path, frame_id: str) -> Optional[Path]:
    """Finds the prepared receptor structure for a given frame."""
    candidates = [
        cwd / "01_prep_receptor" / "structures" / f"receptor_prepared_f{frame_id}.mae",
        cwd / "01_prep_receptor" / "raw" / f"receptor_raw_f{frame_id}.mae"
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def _resolve_ligand_file(cwd: Path, frame_id: str, is_native: bool) -> Optional[Path]:
    """Locates the 3D docked pose or crystal reference file for a ligand."""
    if is_native:
        candidates = [
            cwd / "04a_redocking" / "poses" / f"glide_dock_f{frame_id}_pv.maegz",
            cwd / "01_prep_receptor" / "structures" / f"ref_native_f{frame_id}.mae",
            cwd / "04a_redocking" / "export" / f"glide_dock_f{frame_id}_poses.sdf",
            cwd / "04a_redocking" / "poses" / f"gnina_dock_f{frame_id}.sdf"
        ]
    else:
        candidates = [
            cwd / "04b_screening" / "poses" / f"glide_dock_f{frame_id}_pv.maegz",
            cwd / "04b_screening" / "export" / f"glide_dock_f{frame_id}_poses.sdf",
            cwd / "04b_screening" / "export" / f"gnina_minimize_f{frame_id}_poses.sdf",
            cwd / "04b_screening" / "export" / f"gnina_score_only_f{frame_id}_poses.sdf",
            cwd / "04b_screening" / "export" / f"gnina_dock_f{frame_id}_poses.sdf",
            cwd / "04b_screening" / "poses" / f"gnina_dock_f{frame_id}.sdf"
        ]

    for c in candidates:
        if c.exists():
            return c
    return None


def generate_lid_diagrams(df_resolved: pd.DataFrame, cwd: Path, out_dir: Path, top_n: int = 20, force: bool = False):
    """
    Builds task manifest and executes Schrodinger worker to generate 2D LIDs.
    Selects Top-N library hits, all benchmark controls, and native reference.
    """
    schrodinger_path = os.environ.get("SCHRODINGER")
    if not schrodinger_path:
        logger.warning("SCHRODINGER environment variable is not set. Skipping LID generation.")
        return

    fdir_env = os.environ.get("SBDD_FRAMEWORK_DIR")
    if fdir_env:
        worker_script = Path(fdir_env) / "bin" / "generate_lids.py"
        if not worker_script.exists():
            worker_script = Path(fdir_env) / "src" / "schrodinger_workers" / "generate_lids.py"
    else:
        worker_script = Path(__file__).resolve().parents[1] / "generate_lids.py"

    if not worker_script.exists():
        logger.error(f"LID worker script not found at {worker_script}")
        return

    lids_dir = out_dir / "lids"
    lids_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = lids_dir / "lid_manifest.json"

    # Select candidates: Native + Controls + Top-N Library Hits
    native_df = df_resolved[df_resolved["Is_Native"] == True]
    bench_df = df_resolved[df_resolved.get("Is_Benchmark", False) == True]
    library_top = df_resolved[(df_resolved["Is_Native"] == False) & (~df_resolved["Title"].isin(bench_df["Title"]))].head(top_n)

    target_df = pd.concat([native_df, bench_df, library_top], ignore_index=True).drop_duplicates(subset="Title")

    if target_df.empty:
        logger.warning("No targets available for 2D LID generation.")
        return

    manifest_tasks: List[Dict[str, Any]] = []

    # Map rank among library hits for standardized naming
    lib_titles_ranked = df_resolved[(df_resolved["Is_Native"] == False) & (~df_resolved["Title"].isin(bench_df["Title"]))]["Title"].tolist()

    for _, row in target_df.iterrows():
        title = str(row["Title"])
        frame_id = str(row["Frame"]).replace("F", "")
        is_native = bool(row.get("Is_Native", False))
        is_bench = bool(row.get("Is_Benchmark", False))

        rec_file = _resolve_receptor_file(cwd, frame_id)
        lig_file = _resolve_ligand_file(cwd, frame_id, is_native)

        if not rec_file or not lig_file:
            logger.warning(f"Could not resolve 3D structure files for '{title}' (Frame {frame_id}). Skipping LID.")
            continue

        safe_title = "".join(c if c.isalnum() or c in "-_" else "_" for c in title)

        if is_native:
            png_name = f"native_ref_f{frame_id}_2d_lid.png"
        elif is_bench:
            png_name = f"control_{safe_title}_f{frame_id}_2d_lid.png"
        else:
            rank_num = lib_titles_ranked.index(title) + 1 if title in lib_titles_ranked else 99
            png_name = f"Rank_{rank_num:02d}_{safe_title}_f{frame_id}_2d_lid.png"

        out_png = lids_dir / png_name

        manifest_tasks.append({
            "title": title,
            "frame": frame_id,
            "is_native": is_native,
            "receptor_file": str(rec_file.resolve()),
            "ligand_file": str(lig_file.resolve()),
            "out_png": str(out_png.resolve())
        })

    if not manifest_tasks:
        logger.error("No valid LID tasks generated.")
        return

    # Checkpoint validation
    all_pngs_exist = all(Path(t["out_png"]).exists() for t in manifest_tasks)
    if all_pngs_exist and not force:
        logger.info(f"All {len(manifest_tasks)} 2D LID diagrams already exist in {lids_dir.name}/. Skipping. Use --force to override.")
        return

    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest_tasks, f, indent=2)

    logger.info(f"Generating 2D Ligand Interaction Diagrams for {len(manifest_tasks)} compound(s)...")

    cmd = [
        f"{schrodinger_path}/run", "python3", str(worker_script),
        "--manifest", str(manifest_path.resolve()),
        "--out_dir", str(lids_dir.resolve())
    ]

    success = run_command(cmd, out_dir, "05_generate_lids")

    if success:
        generated_count = len(list(lids_dir.glob("*.png")))
        logger.success(f"2D Interaction Diagrams complete! Saved {generated_count} image(s) to: {lids_dir.name}/")
    else:
        logger.error("LID worker execution encountered errors. Inspect logs/05_generate_lids.log.")