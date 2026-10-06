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


def add_pocket_group(ndx_path: Path, pdb_path: Path, ligand_name: str):
    """Selects protein residues within 5.0 A of ligand plus ligand itself."""
    import MDAnalysis as mda

    u = mda.Universe(str(pdb_path))
    lig_res = (ligand_name or "KVR").strip().upper()
    lig_atoms = u.select_atoms(f"resname {lig_res}")
    if len(lig_atoms) == 0:
        pocket = u.select_atoms("protein")
    else:
        pocket = u.select_atoms(
            f"(same residue as (protein and around 5.0 resname {lig_res})) or resname {lig_res}"
        )
        if len(pocket) == 0:
            pocket = u.select_atoms("protein")

    atom_ids = list(map(int, pocket.atoms.ids))
    with open(ndx_path, "a") as f:
        write_ndx_group(f, "Pocket", atom_ids)
    logging.info(
        f"Injected [ Pocket ] group ({len(atom_ids)} atoms around {lig_res}) into {ndx_path}"
    )


def add_thermostat_groups_by_name(ndx_path: Path, water_names: str, ion_names: str):
    """
    Builds PROT_LIG and SOL_ION groups strictly by string matching config names.
    """
    groups = parse_ndx(ndx_path)

    if "System" not in groups:
        raise KeyError("Group [ System ] not found in index.ndx")

    water_and_ion_atoms: Set[int] = set()
    w_set = {x.strip() for x in water_names.split(",") if x.strip()}
    i_set = {x.strip() for x in ion_names.split(",") if x.strip()}

    matched_water = False
    for w_name in w_set:
        if w_name in groups:
            water_and_ion_atoms.update(groups[w_name])
            matched_water = True
            
    if not matched_water:
        raise KeyError(f"None of the config water names {w_set} found in index.ndx")

    for ion_name in i_set:
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


def generate_posre(top_file: Path, forces_str: str, water_names: str, ion_names: str):
    """Parses topology and dynamically generates position restraints strictly based on user configuration."""
    if not forces_str or not forces_str.strip():
        raise ValueError(
            "FATAL: --forces argument cannot be empty. Restraint forces must be specified explicitly in md.yaml."
        )

    raw_forces = [int(float(x.strip())) for x in forces_str.split(",") if x.strip()]
    if not raw_forces:
        raise ValueError("FATAL: No valid force constants parsed from configuration.")

    sorted_forces = sorted(list(set(raw_forces)), reverse=True)
    max_fc = sorted_forces[0]

    force_levels = []
    for fc in sorted_forces:
        suffix = f"_{fc}"
        force_levels.append((suffix, fc))
        if fc == max_fc:
            force_levels.append(("", fc))

    lines = top_file.read_text().splitlines()
    current_mol = None
    in_atoms = False
    mols = {}

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
                if mass > 1.5:
                    mols[current_mol].append(atom_nr)

    ions = {x.strip() for x in ion_names.split(",") if x.strip()}
    water = {x.strip() for x in water_names.split(",") if x.strip()}

    out_lines = []
    state = "NONE"
    current_mol = None

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("[ moleculetype ]") or stripped.startswith("[ system ]"):
            if current_mol and current_mol not in ions and current_mol not in water:
                is_protein = current_mol.startswith("system") or current_mol.startswith(
                    "Protein"
                )
                prefix = "POSRES" if is_protein else "POSRES_LIG"

                for suffix, fc in force_levels:
                    macro = f"{prefix}{suffix}"
                    out_lines.append(
                        f'#ifdef {macro}\n#include "posre_{current_mol}_{fc}.itp"\n#endif'
                    )

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

        for fc in sorted_forces:
            itp_lines = [
                "[ position_restraints ]",
                "; atom  type      fx      fy      fz",
            ]
            for nr in heavy_atoms:
                itp_lines.append(f"{nr:>6}     1  {fc:6d}  {fc:6d}  {fc:6d}")

            itp_file = Path(f"posre_{mol}_{fc}.itp")
            itp_file.write_text("\n".join(itp_lines) + "\n")

        print(
            f"Dynamically generated restraints for {mol} ({len(heavy_atoms)} atoms) with forces: {sorted_forces}"
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--ndx", required=True)
    parser.add_argument("--top", required=False)
    parser.add_argument("--pdb", required=False)
    parser.add_argument("--ligand", default="KVR")
    parser.add_argument(
        "--action",
        choices=["add_groups", "generate_posre", "ensure_complex", "add_pocket"],
        required=True,
    )
    parser.add_argument("--orgs", default="")
    parser.add_argument(
        "--forces", default="", help="Comma-separated list of force constants from YAML"
    )
    parser.add_argument("--water_names", default="")
    parser.add_argument("--ion_names", default="")
    parser.add_argument("--lipid_names", default="")
    args = parser.parse_args()

    if args.action == "add_groups":
        add_thermostat_groups_by_name(Path(args.ndx), args.water_names, args.ion_names)
    elif args.action == "generate_posre":
        if not args.top:
            raise ValueError("FATAL: --top is required for generate_posre")
        generate_posre(Path(args.top), args.forces, args.water_names, args.ion_names)
    elif args.action == "ensure_complex":
        ensure_complex_ndx_group(Path(args.ndx), args.water_names, args.ion_names, args.lipid_names)
    elif args.action == "add_pocket":
        add_pocket_group(Path(args.ndx), Path(args.pdb), args.ligand)
