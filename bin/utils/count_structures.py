"""
Dynamically counts the number of structures in MAE/MAEGZ/SDF using Schrodinger API.
"""
import os
import subprocess
from pathlib import Path
from loguru import logger

def get_structure_count(file_path: Path) -> int:
    schrodinger_path = os.environ.get("SCHRODINGER")
    cmd = [
        f"{schrodinger_path}/run", "python3", "-c",
        f"from schrodinger.structure import count_structures; print(count_structures(r'{file_path.resolve()}'))"
    ]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        return int(res.stdout.strip())
    except Exception as e:
        logger.warning(f"Could not automatically count structures in {file_path.name}: {e}")
        return 0