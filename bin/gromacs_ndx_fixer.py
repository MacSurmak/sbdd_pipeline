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

    if "System" not in groups:
        raise KeyError("Group [ System ] not found in index.ndx")

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


def generate_posre(top_file: Path):
    """Parses topology and generates strictly local position restraints for heavy atoms."""
    lines = top_file.read_text().splitlines()
    current_mol = None
    in_atoms = False

    mols = {}  # name -> list of heavy atom local indices

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith("["):
            section = stripped.strip("[]").strip()
            if section == "moleculetype":
                current_mol = "EXPECT_NAME"
                in_atoms = False
            elif section == "atoms":
                in_atoms = True
            else:
                in_atoms = False
            continue

        if current_mol == "EXPECT_NAME":
            current_mol = stripped.split()[0]
            mols[current_mol] = []
            continue

        if in_atoms and current_mol:
            parts = stripped.split()
            if len(parts) >= 8:
                atom_nr = parts[0]
                mass = float(parts[7])
                if mass > 1.5:  # Heavy atom threshold
                    mols[current_mol].append(atom_nr)

    ions = {
        "K+",
        "NA+",
        "CL-",
        "Cl-",
        "Na+",
        "ZN",
        "MG",
        "CA",
        "ZN2",
        "MG2",
        "CA2",
        "K",
        "CL",
        "NA",
    }
    water = {"SOL", "Water", "HOH", "WAT"}

    out_lines = []
    state = "NONE"
    current_mol = None

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("[ moleculetype ]") or stripped.startswith("[ system ]"):
            if current_mol and current_mol not in ions and current_mol not in water:
                out_lines.append(f"#ifdef POSRES")
                out_lines.append(f'#include "posre_{current_mol}.itp"')
                out_lines.append(f"#endif")

            state = "MOL_START" if stripped.startswith("[ moleculetype ]") else "SYSTEM"
            current_mol = None

        elif state == "MOL_START" and stripped and not stripped.startswith(";"):
            current_mol = stripped.split()[0]
            state = "INSIDE_MOL"

        out_lines.append(line)

    top_file.write_text("\n".join(out_lines) + "\n")

    for mol, heavy_atoms in mols.items():
        if mol in ions or mol in water:
            continue

        # Determine macro name
        if mol.startswith("system") or mol.startswith("Protein"):
            macro = "POSRES_PROT_FC"
        else:
            macro = f"POSRES_{mol}_FC"

        itp_lines = ["[ position_restraints ]", "; atom  type      fx      fy      fz"]
        for nr in heavy_atoms:
            itp_lines.append(f"{nr:>6}     1  {macro}  {macro}  {macro}")

        itp_file = Path(f"posre_{mol}.itp")
        itp_file.write_text("\n".join(itp_lines) + "\n")
        print(
            f"Generated {itp_file} with macro {macro} ({len(heavy_atoms)} heavy atoms)"
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--ndx", required=True)
    parser.add_argument("--top", required=False)
    parser.add_argument(
        "--action",
        choices=["add_groups", "generate_posre", "ensure_complex"],
        required=True,
    )
    parser.add_argument("--orgs", default="")
    args = parser.parse_args()

    if args.action == "add_groups":
        add_thermostat_groups_by_name(Path(args.ndx))
    elif args.action == "generate_posre":
        if args.top:
            generate_posre(Path(args.top))
    elif args.action == "ensure_complex":
        ensure_complex_ndx_group(Path(args.ndx))
