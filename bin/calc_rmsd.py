"""
bin/calc_rmsd.py
Calculates exact symmetry-aware, in-situ heavy-atom RMSD between reference ligand
and docked poses.
Single deterministic implementation via topological graph isomorphism and symmetry
automorphism permutation. Strictly invariant to atom ordering, file format differences,
and bond-type assignments.
"""

import argparse
import csv
import logging
import math
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def read_heavy_atoms_mol(path: Path) -> list:
    """Reads 3D structures from SDF/PDB and returns list of heavy-atom molecules."""
    raw_mols = []
    if path.suffix.lower() == ".sdf":
        # Strict reading: sanitization MUST be enabled to ensure proper valences for template assignment
        suppl = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=True)
        for idx, m in enumerate(suppl, start=1):
            if m is None:
                raise ValueError(
                    f"FATAL: Molecule at index {idx} in {path.name} failed RDKit sanitization (broken valence/coordinates)."
                )
            raw_mols.append(m)
    elif path.suffix.lower() == ".pdb":
        m = Chem.MolFromPDBFile(str(path), sanitize=True, removeHs=True)
        if m is None:
            raise ValueError(
                f"FATAL: Failed to parse and sanitize PDB molecule from {path.name}"
            )
        raw_mols.append(m)
    else:
        raise ValueError(f"Unsupported molecule format: {path.suffix}")

    clean_mols = []
    for m in raw_mols:
        # Ignore protein or solvent residue artifacts
        if m.GetNumAtoms() > 200:
            continue
        if m.GetNumAtoms() > 0:
            clean_mols.append(m)

    return clean_mols


def calculate_pair_rmsd(m_ref, m_dock) -> float:
    """
    Computes exact in-situ Cartesian RMSD across all symmetry automorphisms.
    Assigns reference bond orders to the docked pose to mathematically bypass
    docking engine valence/aromaticity corruption, eliminating combinatoric explosion.
    """
    n_atoms = m_ref.GetNumAtoms()
    if n_atoms != m_dock.GetNumAtoms():
        raise ValueError(
            f"Atom count mismatch: reference has {n_atoms} heavy atoms, docked pose has {m_dock.GetNumAtoms()}"
        )

    try:
        m_dock_fixed = AllChem.AssignBondOrdersFromTemplate(m_ref, m_dock)
    except Exception as e:
        raise ValueError(
            f"Topological mismatch: unable to assign reference bonds to dock pose. {e}"
        )

    matches = m_dock_fixed.GetSubstructMatches(m_ref, uniquify=False)
    if not matches:
        raise ValueError("Failed to obtain substructure isomorphism matches.")

    best_sq_dist = float("inf")
    conf_ref = m_ref.GetConformer()
    conf_dock = m_dock_fixed.GetConformer()

    for match in matches:
        sq_dist = 0.0
        for i in range(n_atoms):
            p_ref = conf_ref.GetAtomPosition(i)
            p_dock = conf_dock.GetAtomPosition(match[i])
            dx = p_ref.x - p_dock.x
            dy = p_ref.y - p_dock.y
            dz = p_ref.z - p_dock.z
            sq_dist += dx * dx + dy * dy + dz * dz

        if sq_dist < best_sq_dist:
            best_sq_dist = sq_dist

    return math.sqrt(best_sq_dist / n_atoms)


def main():
    parser = argparse.ArgumentParser(description="Calculate In-Situ Heavy-Atom RMSD")
    parser.add_argument(
        "--reference", required=True, help="Path to reference ligand (.sdf or .pdb)"
    )
    parser.add_argument(
        "--poses", required=True, help="Path to docked ligand poses (.sdf)"
    )
    parser.add_argument("--receptor", default="", help="Receptor name")
    parser.add_argument("--title", default="", help="Ligand title")
    parser.add_argument(
        "--engine",
        default="generic",
        choices=["glide", "gnina", "generic"],
        help="Docking engine name",
    )
    parser.add_argument("--output", required=True, help="Path to output CSV")
    args = parser.parse_args()

    ref_path = Path(args.reference).resolve()
    poses_path = Path(args.poses).resolve()
    out_path = Path(args.output).resolve()

    if not ref_path.exists():
        raise FileNotFoundError(f"Reference file does not exist: {ref_path}")
    if not poses_path.exists():
        raise FileNotFoundError(f"Docked poses file does not exist: {poses_path}")

    ref_mols = read_heavy_atoms_mol(ref_path)
    if not ref_mols:
        raise ValueError(f"Could not read reference molecule from {ref_path}")
    m_ref = ref_mols[0]

    docked_mols = read_heavy_atoms_mol(poses_path)
    if not docked_mols:
        raise ValueError(f"Could not read docked poses from {poses_path}")

    matched_title = args.title
    best_rmsd = float("inf")

    # Evaluate all docked poses and find the minimum in-situ RMSD
    for idx, m_dock in enumerate(docked_mols, start=1):
        if not matched_title and m_dock.HasProp("_Name"):
            matched_title = m_dock.GetProp("_Name").strip()

        # Strict Fail-Fast: Any topological deviation or missing atom MUST crash the process
        rmsd_val = calculate_pair_rmsd(m_ref, m_dock)
        if rmsd_val < best_rmsd:
            best_rmsd = rmsd_val

    if math.isinf(best_rmsd):
        raise ValueError(
            f"FATAL: No poses could be evaluated in {poses_path.name} against {ref_path.name}"
        )

    val = round(best_rmsd, 3)
    engine_col = (
        "InSitu_RMSD_Glide_A"
        if args.engine == "glide"
        else ("InSitu_RMSD_GNINA_A" if args.engine == "gnina" else "InSitu_RMSD_A")
    )

    final_title = matched_title or args.title
    if not final_title:
        raise ValueError(
            "FATAL: Ligand title could not be resolved from file properties or CLI arguments."
        )

    row = {
        "Receptor": args.receptor,
        "Title": final_title,
        engine_col: val,
        "InSitu_RMSD_A": val,
        "RMSD_Heavy": val,
    }

    fieldnames = ["Receptor", "Title", engine_col, "InSitu_RMSD_A", "RMSD_Heavy"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)

    logging.info(f"Calculated in-situ RMSD ({args.engine}): {val} Å -> {out_path}")


if __name__ == "__main__":
    main()
