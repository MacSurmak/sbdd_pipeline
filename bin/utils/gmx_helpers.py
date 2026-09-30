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
    """
    Builds [ COMPLEX ] strictly containing Protein + Ligands + all Cofactors.
    Excludes water (SOL), bulk ions (NA, CL, K, etc.), and membrane lipids.
    """
    if not ndx_path.exists():
        raise FileNotFoundError(f"FATAL: index.ndx not found at {ndx_path}")

    groups = parse_ndx(ndx_path)
    if "System" not in groups:
        raise KeyError("Group [ System ] not found in index.ndx")

    solvent_and_lipid_names = {
        "SOL",
        "Water",
        "HOH",
        "WAT",
        "TIP3",
        "OPC",
        "NA",
        "CL",
        "K",
        "POT",
        "SOD",
        "CLA",
        "CAL",
        "MG",
        "ZN",
        "CA",
        "NA+",
        "CL-",
        "K+",
        "MG2+",
        "CA2+",
        "ZN2+",
        "POPC",
        "POPE",
        "CHL1",
        "POPS",
        "DPSM",
        "PAPC",
    }

    solvent_atoms: set[int] = set()
    for name in solvent_and_lipid_names:
        if name in groups:
            solvent_atoms.update(groups[name])

    system_atoms = set(groups["System"])
    complex_atoms = sorted(list(system_atoms - solvent_atoms))

    if not complex_atoms:
        # Fallback to Protein if calculation yielded empty set
        complex_atoms = groups.get("Protein", [])

    # Replace or append [ COMPLEX ]
    lines = ndx_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    out_lines = []
    skip = False
    for line in lines:
        if line.strip().lower() in ["[ complex ]"]:
            skip = True
            continue
        if skip and line.strip().startswith("["):
            skip = False
        if not skip:
            out_lines.append(line)

    ndx_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    with open(ndx_path, "a", encoding="utf-8") as f:
        write_ndx_group(f, "COMPLEX", complex_atoms)

    logging.info(
        f"Injected complete [ COMPLEX ] group ({len(complex_atoms)} atoms: Protein + Ligand + Cofactors) into index.ndx"
    )
