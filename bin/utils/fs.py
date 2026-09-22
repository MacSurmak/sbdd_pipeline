"""
File: src/utils/fs.py
Description: File system utilities, including robust numerical sorting for frames.
"""
import re
from pathlib import Path
from typing import List, Tuple

def get_frames_sorted(directory: Path, pattern: str) -> List[Tuple[Path, str]]:
    """
    Finds files by pattern and sorts them numerically by frame ID '_f(\\d+)'.
    Returns a list of tuples: (Path, frame_id_as_string)
    """
    if not directory.exists():
        return []
    
    files = list(directory.glob(pattern))
    
    def extract_frame(p: Path) -> int:
        match = re.search(r'_f(\d+)', p.stem)
        return int(match.group(1)) if match else 0

    sorted_files = sorted(files, key=extract_frame)
    
    return [(f, str(extract_frame(f))) for f in sorted_files if extract_frame(f) > 0]