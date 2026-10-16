"""
File: cli/core/chem.py
Description: Chemoinformatics helpers for ligand standardization and identity verification.
Computes Standard InChIKey and Canonical SMILES via RDKit for topological deduplication.
"""

from pathlib import Path
from typing import Tuple, Optional
import sys


def get_mol_identity(source: str | Path) -> Tuple[str, str, str, int]:
    """
    Parses a chemical structure (SDF file, molblock, SMILES or PDB),
    standardizes it, and returns (inchikey, canonical_smiles, mol_name).
    Raises ValueError on invalid chemistry.
    """
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError(
            "RDKit is required for molecular identity and InChIKey calculation."
        )

    mol = None
    mol_name = ""

    if isinstance(source, Path) or (isinstance(source, str) and Path(source).exists()):
        p = Path(source)
        if p.suffix.lower() in [".sdf", ".mol"]:
            suppl = Chem.SDMolSupplier(str(p), removeHs=False)
            mol = next((m for m in suppl if m is not None), None)
            if mol and mol.HasProp("_Name"):
                mol_name = mol.GetProp("_Name").strip()
        elif p.suffix.lower() == ".pdb":
            mol = Chem.MolFromPDBFile(str(p), removeHs=False)
        else:
            raise ValueError(f"Unsupported molecule file format: {p.suffix}")
    elif isinstance(source, str):
        source_clean = source.strip()
        if "\n" in source_clean or "$$$$" in source_clean:
            mol = Chem.MolFromMolBlock(source_clean, removeHs=False)
            if mol and mol.HasProp("_Name"):
                mol_name = mol.GetProp("_Name").strip()
        else:
            mol = Chem.MolFromSmiles(source_clean)
            mol_name = source_clean

    if mol is None:
        raise ValueError(
            "Could not parse valid chemical structure from provided input."
        )

    from rdkit.Chem.MolStandardize import rdMolStandardize

    formal_charge = Chem.GetFormalCharge(mol)

    # Check if entity is a monoatomic ion to preserve its formal oxidation state
    is_monoatomic_ion = mol.GetNumHeavyAtoms() == 1 and mol.GetNumAtoms() == 1

    if is_monoatomic_ion:
        parent_inchikey = Chem.MolToInchiKey(mol)
        parent_smiles = Chem.MolToSmiles(mol, canonical=True)
    else:
        # Neutralize organic parent to establish invariant chemical substance identity
        # across protonation microstates, while preserving individual pose physics
        uncharger = rdMolStandardize.Uncharger()
        heavy_mol = Chem.RemoveHs(mol)
        parent_mol = uncharger.uncharge(heavy_mol)
        parent_inchikey = Chem.MolToInchiKey(parent_mol)
        parent_smiles = Chem.MolToSmiles(parent_mol, canonical=True)

    if not parent_inchikey:
        raise ValueError("Failed to generate InChIKey. Invalid chemical structure.")

    if not mol_name:
        mol_name = parent_inchikey[:14]

    return parent_inchikey, parent_smiles, mol_name, formal_charge
