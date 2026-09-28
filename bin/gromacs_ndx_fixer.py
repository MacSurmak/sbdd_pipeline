#!/usr/bin/env python3
"""
File: src/utils/gromacs_utils.py
Description: Safe topology and index manipulations for GROMACS.
Builds thermostat groups strictly by dictionary name lookup, NO index offsets.
"""
import os
import sys
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Set

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def get_gmx_binary(config: dict = None) -> str:
    """
    Strictly validates and returns the system GROMACS binary 'gmx' from PATH.
    """
    if not shutil.which("gmx"):
        logging.critical("System GROMACS binary 'gmx' not found in PATH.")
        sys.exit(1)
    return "gmx"


# Add script directory to sys.path for absolute imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.gmx_helpers import parse_ndx, write_ndx_group, ensure_complex_ndx_group


def add_thermostat_groups_by_name(ndx_path: Path):
    """
    Builds PROT_LIG and SOL_ION groups strictly by string matching names in index.ndx.
    NO numerical index offsets. Guaranteed partition of the entire System.
    """
    groups = parse_ndx(ndx_path)
    prot_atoms = groups.get("Protein", [])
    unk_atoms = groups.get("UNK", [])
    if not prot_atoms:
        raise ValueError("FATAL: 'Protein' group not found in index.ndx")

    complex_atoms = sorted(list(set(prot_atoms + unk_atoms)))

    water_and_ion_atoms: Set[int] = set()
    ion_names = {
        "K+",
        "Cl-",
        "NA+",
        "CL-",
        "NA",
        "CL",
        "K",
        "POT",
        "CLA",
        "SOD",
        "MG",
        "ZN",
        "CA",
        "CAL",
    }

    # 1. Match water strictly by name
    if "SOL" in groups:
        water_and_ion_atoms.update(groups["SOL"])
    elif "Water" in groups:
        water_and_ion_atoms.update(groups["Water"])
    else:
        raise KeyError("Neither [ SOL ] nor [ Water ] found in index.ndx")

    # 2. Match ions strictly by name
    for ion_name in ion_names:
        if ion_name in groups:
            water_and_ion_atoms.update(groups[ion_name])

    # 3. Everything else in System is solute + membrane (Protein + UNK + Lipids)
    system_atoms = set(groups["System"])
    solute_membrane_atoms = system_atoms - water_and_ion_atoms

    # Strict partition validation
    if not solute_membrane_atoms:
        raise ValueError("PROT_LIG group is empty! Check index.ndx names.")
    if not water_and_ion_atoms:
        raise ValueError("SOL_ION group is empty! Check index.ndx names.")
    if len(solute_membrane_atoms) + len(water_and_ion_atoms) != len(system_atoms):
        raise ValueError("Partition mismatch: PROT_LIG + SOL_ION != System!")

    # Append new groups strictly formatted
    with open(ndx_path, "a") as f:
        write_ndx_group(f, "PROT_LIG", sorted(list(solute_membrane_atoms)))
        write_ndx_group(f, "SOL_ION", sorted(list(water_and_ion_atoms)))


def inject_ligand_posre(top_file: Path):
    """Injects position restraints with multiple force constants at the end of UNK moleculetype."""
    if not top_file.exists():
        return
    lines = top_file.read_text().splitlines()

    out_lines = []
    in_unk = False
    injected = False

    posre_block = (
        "\n#ifdef POSRES_LIG\n"
        '#include "posre_lig.itp"\n'
        "#endif\n"
        "#ifdef POSRES_LIG_250\n"
        '#include "posre_lig_250.itp"\n'
        "#endif\n"
        "#ifdef POSRES_LIG_50\n"
        '#include "posre_lig_50.itp"\n'
        "#endif\n"
    )

    for line in lines:
        stripped = line.strip()

        if (
            in_unk
            and not injected
            and (
                stripped.startswith("[ moleculetype ]")
                or stripped.startswith("[ system ]")
            )
        ):
            out_lines.append(posre_block)
            injected = True
            in_unk = False

        if stripped.startswith("[ moleculetype ]"):
            in_unk = False

        if (
            not stripped.startswith(";")
            and stripped.startswith("UNK")
            and len(stripped.split()) >= 2
        ):
            in_unk = True

        out_lines.append(line)

    if in_unk and not injected:
        out_lines.append(posre_block)

    # Inject protein multi-level posres at the end of system1 or Protein
    out_lines_prot = []
    in_prot = False
    injected_prot = False

    posre_block_prot = (
        "\n#ifdef POSRES\n"
        '#include "posre.itp"\n'
        "#endif\n"
        "#ifdef POSRES_250\n"
        '#include "posre_250.itp"\n'
        "#endif\n"
        "#ifdef POSRES_50\n"
        '#include "posre_50.itp"\n'
        "#endif\n"
    )

    for line in out_lines:
        stripped = line.strip()

        if (
            in_prot
            and not injected_prot
            and (
                stripped.startswith("[ moleculetype ]")
                or stripped.startswith("[ system ]")
            )
        ):
            out_lines_prot.append(posre_block_prot)
            injected_prot = True
            in_prot = False

        if stripped.startswith("[ moleculetype ]"):
            in_prot = False

        if (
            not stripped.startswith(";")
            and (stripped.startswith("system1") or stripped.startswith("Protein"))
            and len(stripped.split()) >= 2
        ):
            in_prot = True

        out_lines_prot.append(line)

    top_file.write_text("\n".join(out_lines_prot) + "\n")


def fix_posre_indices(itp_file: Path):
    """Converts global atom indices in posre_lig.itp to 1-based local molecule indices (1 to N)."""
    if not itp_file.exists():
        return
    lines = itp_file.read_text().splitlines()
    out_lines = []
    atom_idx = 1
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(";") or stripped.startswith("["):
            out_lines.append(line)
            continue
        parts = stripped.split()
        if len(parts) >= 5:
            out_lines.append(
                f"{atom_idx:5d}{int(parts[1]):5d}{float(parts[2]):10.1f}{float(parts[3]):10.1f}{float(parts[4]):10.1f}"
            )
            atom_idx += 1
        else:
            out_lines.append(line)
    itp_file.write_text("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--ndx", required=True)
    parser.add_argument("--top", required=False)
    parser.add_argument(
        "--action",
        choices=["add_groups", "inject_posre", "fix_posre", "ensure_complex"],
        required=True,
    )
    args = parser.parse_args()

    ndx_path = Path(args.ndx)
    if args.action == "add_groups":
        add_thermostat_groups_by_name(ndx_path)
    elif args.action == "inject_posre":
        if args.top:
            inject_ligand_posre(Path(args.top))
    elif args.action == "fix_posre":
        for sfx in ["", "_250", "_50"]:
            fix_posre_indices(Path(f"posre{sfx}.itp"))
            fix_posre_indices(Path(f"posre_lig{sfx}.itp"))
    elif args.action == "ensure_complex":
        ensure_complex_ndx_group(ndx_path)
