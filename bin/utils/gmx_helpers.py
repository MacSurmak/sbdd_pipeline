"""
Common GROMACS index and topology helper routines.
"""

from pathlib import Path
from typing import Dict, List
import logging


def parse_ndx(ndx_path: Path) -> Dict[str, List[int]]:
    """Parses an NDX file into a dictionary of group_name -> list of atom IDs."""
    groups: Dict[str, List[int]] = {}
    current_group = None
    with open(ndx_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                current_group = line[1:-1].strip()
                groups[current_group] = []
            elif current_group is not None:
                groups[current_group].extend([int(x) for x in line.split()])
    return groups


def write_ndx_group(f, name: str, atom_ids: List[int]):
    """Writes a single NDX group strictly formatted in columns of 15."""
    f.write(f"\n[ {name} ]\n")
    for i in range(0, len(atom_ids), 15):
        f.write(" ".join(f"{x:5d}" for x in atom_ids[i : i + 15]) + "\n")


def ensure_complex_ndx_group(ndx_path: Path):
    """Parses index.ndx and injects [ COMPLEX ] (Protein + UNK) if missing."""
    if not ndx_path.exists():
        return
    text = ndx_path.read_text(encoding="utf-8", errors="ignore")
    if "[ COMPLEX ]" in text or "[ complex ]" in text:
        return

    groups = parse_ndx(ndx_path)
    prot_atoms = groups.get("Protein", [])
    unk_atoms = groups.get("UNK", [])
    if not prot_atoms:
        logging.warning("'Protein' group not found in index.ndx")
        return

    complex_atoms = sorted(list(set(prot_atoms + unk_atoms)))
    with open(ndx_path, "a", encoding="utf-8") as f:
        write_ndx_group(f, "COMPLEX", complex_atoms)
    logging.info(
        f"Injected [ COMPLEX ] group ({len(complex_atoms)} atoms) into index.ndx"
    )
