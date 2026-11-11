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


def count_dry_top_atoms(dry_top_path: Path) -> int:
    """
    Deterministically counts the exact number of atoms in the dry complex topology
    (complex_dry.top or complex_dry.prmtop) to preserve structural pocket ions.
    """
    if not dry_top_path.exists() or dry_top_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Dry topology file missing or empty: {dry_top_path}"
        )
    import parmed as pmd

    struct = pmd.load_file(str(dry_top_path))
    if len(struct.atoms) == 0:
        raise ValueError(f"FATAL: Dry topology {dry_top_path} resolved to 0 atoms!")
    return len(struct.atoms)


def ensure_complex_ndx_group(
    ndx_path: Path,
    water_names: str,
    ion_names: str,
    lipid_names: str,
    dry_top_path: Path | None = None,
):
    """
    Builds [ COMPLEX ] strictly containing Protein + Ligands + all structural Cofactors/Ions.
    When dry_top_path is provided, guarantees exact atom count parity with complex_dry.top
    so structural ions in the binding pocket are never stripped along with bulk salt.
    """
    if not ndx_path.exists():
        raise FileNotFoundError(f"FATAL: index.ndx not found at {ndx_path}")

    groups = parse_ndx(ndx_path)
    if "System" not in groups:
        raise KeyError("Group [ System ] not found in index.ndx")

    if dry_top_path is not None:
        n_dry = count_dry_top_atoms(dry_top_path)
        system_atoms_sorted = sorted(groups["System"])
        if len(system_atoms_sorted) < n_dry:
            raise ValueError(
                f"FATAL: System atom count ({len(system_atoms_sorted)}) is smaller than dry topology ({n_dry})!"
            )
        complex_atoms = system_atoms_sorted[:n_dry]
    else:
        w_set = {x.strip() for x in water_names.split(",") if x.strip()}
        i_set = {x.strip() for x in ion_names.split(",") if x.strip()}
        l_set = {x.strip() for x in lipid_names.split(",") if x.strip()}
        solvent_and_lipid_names = w_set | i_set | l_set

        solvent_atoms: set[int] = set()
        for name in solvent_and_lipid_names:
            if name in groups:
                solvent_atoms.update(groups[name])

        system_atoms = set(groups["System"])
        complex_atoms = sorted(list(system_atoms - solvent_atoms))

    if not complex_atoms:
        raise ValueError(
            "FATAL: Group [ COMPLEX ] is empty after excluding solvent and lipids! "
            "Check water_names, ion_names, and lipid_names in md.yaml configuration."
        )

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
