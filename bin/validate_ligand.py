#!/usr/bin/env python3
"""
bin/validate_ligand.py
Fast pre-flight validation of molecular 3D geometry and spin multiplicity calculation.
Strict Fail-Fast: rejects degenerate coordinates, unphysical bonds, and missing hydrogens.
"""

import argparse
import json
import math
import sys
from pathlib import Path
from rdkit import Chem


def validate_geometry(mol: Chem.Mol, mol_id: str, source_path: Path):
    """Enforces physical sanity checks on 3D coordinates."""
    has_hydrogens = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
    if not has_hydrogens:
        raise ValueError(
            f"FATAL: Molecule '{mol_id}' in {source_path} has no explicit hydrogens! "
            "Quantum parametrization requires a fully protonated 3D structure."
        )

    conf = mol.GetConformer()
    for idx, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(idx)
        if abs(pos.x) < 1e-4 and abs(pos.y) < 1e-4 and abs(pos.z) < 1e-4:
            raise ValueError(
                f"FATAL: Atom {atom.GetSymbol()}{idx+1} in molecule '{mol_id}' has degenerate "
                f"coordinates (0, 0, 0) in {source_path}."
            )

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        p1 = conf.GetAtomPosition(a1.GetIdx())
        p2 = conf.GetAtomPosition(a2.GetIdx())
        dist = math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)
        if dist < 0.6:
            raise ValueError(
                f"FATAL: Unphysical steric clash between {a1.GetSymbol()}{a1.GetIdx()+1} and "
                f"{a2.GetSymbol()}{a2.GetIdx()+1}: {dist:.3f} Å < 0.6 Å in {source_path}"
            )
        if (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1) and dist > 1.8:
            raise ValueError(
                f"FATAL: Severed hydrogen bond between {a1.GetSymbol()}{a1.GetIdx()+1} and "
                f"{a2.GetSymbol()}{a2.GetIdx()+1}: {dist:.3f} Å > 1.8 Å in {source_path}"
            )


def calculate_spin_multiplicity(mol: Chem.Mol, formal_charge: int) -> int:
    """Calculates deterministic spin multiplicity (2S + 1) based on valence electrons."""
    total_protons = sum(atom.GetAtomicNum() for atom in mol.GetAtoms())
    total_electrons = total_protons - formal_charge
    return 2 if (total_electrons % 2 != 0) else 1


def main():
    parser = argparse.ArgumentParser(
        description="Validate molecular 3D geometry and determine spin multiplicity"
    )
    parser.add_argument("--sdf", required=True, help="Input SDF structure")
    parser.add_argument("--mol_id", required=True, help="Molecule identifier")
    parser.add_argument("--charge", required=True, type=int, help="Net formal charge")
    parser.add_argument("--out_sdf", required=True, help="Validated output SDF")
    parser.add_argument(
        "--out_json", required=True, help="Output metadata JSON deliverable"
    )
    args = parser.parse_args()

    sdf_path = Path(args.sdf).resolve()
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Input SDF deliverable not found: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(
            f"FATAL: Could not read valid molecular structure from {sdf_path}"
        )

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        mol.UpdatePropertyCache(strict=False)

    validate_geometry(mol, args.mol_id, sdf_path)
    mult = calculate_spin_multiplicity(mol, args.charge)

    mol.SetProp("SpinMultiplicity", str(mult))
    mol.SetProp("FormalCharge", str(args.charge))
    mol.SetProp("_Name", args.mol_id)

    out_sdf_path = Path(args.out_sdf).resolve()
    out_sdf_path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(out_sdf_path)) as writer:
        writer.SetForceV3000(False)
        writer.write(mol)

    metadata = {
        "mol_id": args.mol_id,
        "num_atoms": mol.GetNumAtoms(),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "formal_charge": args.charge,
        "spin_multiplicity": mult,
    }

    out_json_path = Path(args.out_json).resolve()
    out_json_path.parent.mkdir(parents=True, exist_ok=True)
    out_json_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(
        f"✓ Validated ligand geometry for '{args.mol_id}' (Atoms: {mol.GetNumAtoms()}, Multiplicity: {mult})",
        flush=True,
    )


if __name__ == "__main__":
    main()
