"""
File: src/utils/schrodinger.py
Description: Wrappers for Schrodinger Python API workers to avoid code duplication.
"""
import os
import csv
from pathlib import Path
from loguru import logger
from utils.cmd_runner import run_command

def run_rmsd_worker(ref_ligand: Path, poses_file: Path, out_csv: Path, work_dir: Path) -> tuple:
    """
    Invokes the run_rmsd.py worker and parses the best pose result.
    Returns (best_rmsd, best_score) or (None, None).
    """
    schrodinger_path = os.environ.get("SCHRODINGER")
    framework_dir = Path(os.environ.get("SBDD_FRAMEWORK_DIR"))
    worker_script = framework_dir / "src" / "schrodinger_workers" / "run_rmsd.py"
    
    cmd_rmsd = [
        f"{schrodinger_path}/run", "python3", str(worker_script),
        "--reference", str(ref_ligand.resolve()),
        "--poses", str(poses_file.resolve()),
        "--output", str(out_csv.resolve())
    ]
    
    run_command(cmd_rmsd, work_dir, f"rmsd_{poses_file.stem}")
    
    if out_csv.exists():
        with open(out_csv, "r") as f_rmsd:
            rows = list(csv.DictReader(f_rmsd))
        if rows:
            best_pose = min(rows, key=lambda r: float(r["RMSD_Heavy"]))
            score_val = float(best_pose.get("Score", 0.0) or best_pose.get("GlideScore", 0.0) or 0.0)
            rmsd_val = float(best_pose["RMSD_Heavy"])
            return rmsd_val, score_val
            
    return None, None