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

    if hasattr(source, "GetAtoms"):
        mol = Chem.Mol(source)
        if mol.HasProp("_Name"):
            mol_name = mol.GetProp("_Name").strip()
    elif isinstance(source, Path) or (
        isinstance(source, str) and Path(source).exists()
    ):
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


def get_microstate_identity(source: str | Path) -> Tuple[str, str, str, str, str, int]:
    """
    Computes both the Abstract Substance Identity (Neutral InChIKey) and the
    Exact Chemical Microstate Identity (Isomeric Canonical SMILES and state hash).

    Returns:
        (parent_inchikey, parent_smiles, mol_name, state_hash, isomeric_smiles, formal_charge)
    """
    import hashlib

    parent_inchikey, parent_smiles, mol_name, formal_charge = get_mol_identity(source)

    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("RDKit is required for microstate identification.")

    mol = None
    if hasattr(source, "GetAtoms"):
        mol = Chem.Mol(source)
    elif isinstance(source, Path) or (
        isinstance(source, str) and Path(source).exists()
    ):
        p = Path(source)
        if p.suffix.lower() in [".sdf", ".mol"]:
            suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=False)
            mol = next((m for m in suppl if m is not None), None)
        elif p.suffix.lower() == ".pdb":
            mol = Chem.MolFromPDBFile(str(p), removeHs=False, sanitize=False)
    elif isinstance(source, str):
        source_clean = source.strip()
        if "\n" in source_clean or "$$$$" in source_clean:
            mol = Chem.MolFromMolBlock(source_clean, removeHs=False, sanitize=False)
        else:
            mol = Chem.MolFromSmiles(source_clean, sanitize=False)

    if mol is None:
        raise ValueError(
            f"Could not parse valid chemical structure from source: {source}"
        )

    # Standardize property cache to ensure valences and formal charges are evaluated
    mol.UpdatePropertyCache(strict=False)

    # Compute Isomeric Canonical SMILES with stereochemistry and explicit charges
    isomeric_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    state_hash = hashlib.sha256(isomeric_smiles.encode("utf-8")).hexdigest()[:12]

    return (
        parent_inchikey,
        parent_smiles,
        mol_name,
        state_hash,
        isomeric_smiles,
        formal_charge,
    )


def signed_chiral_volume(conf, center_idx: int, nbr_indices: list[int]) -> float:
    """
    Computes signed scalar triple product for a 4-coordinate tetrahedral center:
    V = (v1 - v0) . [(v2 - v0) x (v3 - v0)]
    Enforces strict mathematical invariant for R/S stereochemical parity validation.
    """
    import numpy as np

    if len(nbr_indices) < 3:
        return 0.0
    p0 = conf.GetAtomPosition(center_idx)
    p1 = conf.GetAtomPosition(nbr_indices[0])
    p2 = conf.GetAtomPosition(nbr_indices[1])
    p3 = conf.GetAtomPosition(nbr_indices[2])
    v1 = np.array([p1.x - p0.x, p1.y - p0.y, p1.z - p0.z], dtype=np.float64)
    v2 = np.array([p2.x - p0.x, p2.y - p0.y, p2.z - p0.z], dtype=np.float64)
    v3 = np.array([p3.x - p0.x, p3.y - p0.y, p3.z - p0.z], dtype=np.float64)
    return float(np.dot(v1, np.cross(v2, v3)))


def transfer_coordinates_to_template(
    template_mol,
    source_mol,
    relax_missing_hydrogens: bool = False,
):
    """
    Universal Single Source of Truth for Chemical Topology Restoration & Coordinate Transfer.
    
    Architectural Contract:
    - Coordinates are taken strictly from source_mol (e.g. docked pose, MD frame).
    - Physical topology (bond orders, formal charges, stereochemistry, valences) is taken
      strictly and immutably from template_mol (canonical SubstanceState).
    - Evaluates all topological symmetry automorphisms (handling 180 deg benzene ring flips,
      carboxylate/phosphate oxygen permutations) to select the minimum in-situ heavy-atom RMSD.
    - Validates tetrahedral stereocenter parity via signed chiral volume (Scalar Triple Product).
      Fail-Fast: crashes immediately with fatal ValueError if chirality is inverted.
    - If relax_missing_hydrogens is True (e.g. for GNINA united-atom poses lacking protons),
      heavy atoms are held strictly rigid with infinite force constants (Fixed Points) while
      newly constructed hydrogens are relaxed in the force field (MMFF94 / UFF).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np

    if template_mol is None or source_mol is None:
        raise ValueError("FATAL: transfer_coordinates_to_template received None for template or source molecule!")

    tmpl = Chem.Mol(template_mol)
    tmpl_conf = tmpl.GetConformer()
    src_conf = source_mol.GetConformer()

    tmpl_heavy = Chem.RemoveHs(tmpl, sanitize=False)
    src_heavy = Chem.RemoveHs(source_mol, sanitize=False)

    n_tmpl_heavy = tmpl_heavy.GetNumAtoms()
    n_src_heavy = src_heavy.GetNumAtoms()

    if n_tmpl_heavy != n_src_heavy:
        raise ValueError(
            f"FATAL: Heavy atom count mismatch during coordinate transfer! "
            f"Template has {n_tmpl_heavy} heavy atoms, source has {n_src_heavy}."
        )

    # Element frequency verification
    tmpl_elements = sorted([a.GetAtomicNum() for a in tmpl_heavy.GetAtoms()])
    src_elements = sorted([a.GetAtomicNum() for a in src_heavy.GetAtoms()])
    if tmpl_elements != src_elements:
        raise ValueError(
            f"FATAL: Chemical element composition mismatch between template and source molecule! "
            f"Template: {tmpl_elements} vs Source: {src_elements}"
        )

    # Identify true tetrahedral chiral centers directly on tmpl_heavy to guarantee identical atom indexing
    Chem.AssignStereochemistry(tmpl_heavy, force=True, cleanIt=True)
    chiral_centers = []
    tmpl_heavy_conf = tmpl_heavy.GetConformer()
    for a in tmpl_heavy.GetAtoms():
        if a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            nbrs = [nbr.GetIdx() for nbr in a.GetNeighbors()]
            if len(nbrs) >= 3:
                vol = signed_chiral_volume(tmpl_heavy_conf, a.GetIdx(), nbrs)
                if abs(vol) > 0.05:
                    chiral_centers.append((a.GetIdx(), nbrs, np.sign(vol)))

    # Pure connectivity queries without bond orders to resolve Kekule and thermal bond strain issues
    q_tmpl = Chem.Mol(tmpl_heavy)
    q_src = Chem.Mol(src_heavy)
    for q in (q_tmpl, q_src):
        for b in q.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        for a in q.GetAtoms():
            a.SetFormalCharge(0)
            a.SetIsotope(0)

    matches = q_src.GetSubstructMatches(q_tmpl, uniquify=False)
    if not matches:
        raise ValueError(
            "FATAL: Topological graph isomorphism failed! "
            "Source molecule connectivity does not match template graph."
        )

    src_heavy_conf = src_heavy.GetConformer()
    best_match = None
    best_sq_dist = float("inf")

    # Evaluate all symmetry automorphisms and enforce strict chiral parity
    for match in matches:
        chiral_valid = True
        for tmpl_center, tmpl_nbrs, expected_sign in chiral_centers:
            src_center = match[tmpl_center]
            src_nbrs = [match[n] for n in tmpl_nbrs]
            src_vol = signed_chiral_volume(src_heavy_conf, src_center, src_nbrs)
            if np.sign(src_vol) != expected_sign and abs(src_vol) > 0.05:
                chiral_valid = False
                break

        if not chiral_valid:
            continue

        sq_dist = 0.0
        for tmpl_idx, src_idx in enumerate(match):
            p_t = tmpl_heavy_conf.GetAtomPosition(tmpl_idx)
            p_s = src_heavy_conf.GetAtomPosition(src_idx)
            sq_dist += (p_t.x - p_s.x) ** 2 + (p_t.y - p_s.y) ** 2 + (p_t.z - p_s.z) ** 2

        if sq_dist < best_sq_dist:
            best_sq_dist = sq_dist
            best_match = match

    # Fail-Fast: Zero tolerance for stereochemical inversion
    if best_match is None:
        raise ValueError(
            "FATAL: Stereochemical parity check failed across all topological automorphisms! "
            "Tetrahedral center inversion detected. Execution aborted by Fail-Fast law."
        )

    # Create destination molecule as pristine clone of canonical template
    res_mol = Chem.Mol(tmpl)
    res_conf = res_mol.GetConformer()

    # Map heavy atom coordinates onto result conformer
    for tmpl_heavy_idx, src_heavy_idx in enumerate(best_match):
        orig_tmpl_idx = tmpl_heavy.GetAtomWithIdx(tmpl_heavy_idx).GetIdx()
        orig_src_idx = src_heavy.GetAtomWithIdx(src_heavy_idx).GetIdx()
        pos = src_conf.GetAtomPosition(orig_src_idx)
        if abs(pos.x) < 1e-4 and abs(pos.y) < 1e-4 and abs(pos.z) < 1e-4:
            raise ValueError(
                f"FATAL: Degenerate coordinate (0, 0, 0) detected at atom index {orig_src_idx} during transfer!"
            )
        res_conf.SetAtomPosition(orig_tmpl_idx, pos)

    # Reconstruct/relax missing hydrogens if requested (e.g. for GNINA united-atom poses lacking non-polar protons)
    if relax_missing_hydrogens and source_mol.GetNumAtoms() < res_mol.GetNumAtoms():
        # Re-initialize hydrogen 3D positions based on updated heavy atoms
        temp_h_mol = Chem.RemoveHs(res_mol)
        temp_h_mol = Chem.AddHs(temp_h_mol, addCoords=True)

        mp = AllChem.MMFFGetMoleculeProperties(temp_h_mol)
        ff = AllChem.MMFFGetMoleculeForceField(temp_h_mol, mp) if mp else AllChem.UFFGetMoleculeForceField(temp_h_mol)
        if not ff:
            raise ValueError("FATAL: Could not build force field (MMFF/UFF) to relax missing hydrogens!")

        # Hold all heavy atoms strictly rigid; relax only protons
        for atom in temp_h_mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                ff.AddFixedPoint(atom.GetIdx())
        ff.Minimize(maxIts=500)

        temp_conf = temp_h_mol.GetConformer()
        for atom in res_mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                h_pos = temp_conf.GetAtomPosition(atom.GetIdx())
                res_conf.SetAtomPosition(atom.GetIdx(), h_pos)

    # Validate that no atom in final molecule has degenerate coordinates
    for i in range(res_mol.GetNumAtoms()):
        p = res_conf.GetAtomPosition(i)
        if abs(p.x) < 1e-4 and abs(p.y) < 1e-4 and abs(p.z) < 1e-4:
            sym = res_mol.GetAtomWithIdx(i).GetSymbol()
            raise ValueError(f"FATAL: Transferred molecule contains degenerate coordinate (0, 0, 0) at atom {sym}{i+1}!")

    # Preserve essential metadata properties
    for prop in source_mol.GetPropNames():
        res_mol.SetProp(prop, source_mol.GetProp(prop))
    if source_mol.HasProp("_Name"):
        res_mol.SetProp("_Name", source_mol.GetProp("_Name"))

    return res_mol


def compute_symmetry_rmsd(template_mol, target_mol) -> float:
    """
    Computes exact in-situ heavy-atom RMSD across all symmetry automorphisms.
    Single Source of Truth implementation matching transfer_coordinates_to_template.
    """
    from rdkit import Chem
    import numpy as np

    tmpl_heavy = Chem.RemoveHs(template_mol, sanitize=False)
    target_heavy = Chem.RemoveHs(target_mol, sanitize=False)

    n_atoms = tmpl_heavy.GetNumAtoms()
    if n_atoms != target_heavy.GetNumAtoms():
        raise ValueError(
            f"FATAL: Atom count mismatch for RMSD: template has {n_atoms}, target has {target_heavy.GetNumAtoms()}"
        )

    q_tmpl = Chem.Mol(tmpl_heavy)
    q_target = Chem.Mol(target_heavy)
    for q in (q_tmpl, q_target):
        for b in q.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        for a in q.GetAtoms():
            a.SetFormalCharge(0)
            a.SetIsotope(0)

    matches = q_target.GetSubstructMatches(q_tmpl, uniquify=False)
    if not matches:
        raise ValueError("FATAL: Failed to obtain topological isomorphism matches for RMSD calculation!")

    conf_tmpl = tmpl_heavy.GetConformer()
    conf_target = target_heavy.GetConformer()
    best_sq_dist = float("inf")

    for match in matches:
        sq_dist = 0.0
        for tmpl_idx, target_idx in enumerate(match):
            p_t = conf_tmpl.GetAtomPosition(tmpl_idx)
            p_tgt = conf_target.GetAtomPosition(target_idx)
            sq_dist += (p_t.x - p_tgt.x) ** 2 + (p_t.y - p_tgt.y) ** 2 + (p_t.z - p_tgt.z) ** 2

        if sq_dist < best_sq_dist:
            best_sq_dist = sq_dist

    return float(np.sqrt(best_sq_dist / float(n_atoms)))
