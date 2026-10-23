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
    source_mol_or_coords=None,
    source_mol=None,
    source_elements=None,
    relax_missing_hydrogens: bool = False,
):
    """
    Universal Single Source of Truth for Chemical Topology Restoration & Coordinate Transfer.
    
    Supports two distinct physical input modalities:
    1. NumPy Array (N, 3): Direct in-memory coordinates from MD simulation frames (MDAnalysis).
       Executes 0.01 ms direct index transfer with element verification and built-in physical bond gating.
    2. RDKit Chem.Mol: 3D conformations from docking engines (Glide, GNINA).
       Evaluates all symmetry automorphisms, enforces chiral volume parity (V_chiral), selects minimal
       in-situ RMSD, and reconstructs missing non-polar hydrogens (MMFF94/UFF) with fixed heavy atoms.

    Built-in Fail-Fast physical assertion:
    Verifies that all covalent X-H bond lengths are strictly within 0.80 - 1.35 Å and heavy bonds within 1.00 - 2.20 Å.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np

    src = source_mol_or_coords if source_mol_or_coords is not None else source_mol
    if template_mol is None or src is None:
        raise ValueError("FATAL: transfer_coordinates_to_template received None for template or source!")

    tmpl = Chem.Mol(template_mol)
    tmpl_n_atoms = tmpl.GetNumAtoms()
    res_mol = Chem.Mol(tmpl)
    res_conf = res_mol.GetConformer()

    # -------------------------------------------------------------------------
    # MODALITY 1: Direct NumPy / Sequence Coordinate Array (All-Atom MD Stream)
    # -------------------------------------------------------------------------
    if isinstance(src, (np.ndarray, list)):
        coords = np.asarray(src, dtype=np.float64)
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(f"FATAL: Expected coordinate array of shape (N, 3), got: {coords.shape}")

        if coords.shape[0] != tmpl_n_atoms:
            raise ValueError(
                f"FATAL: Coordinate count mismatch! Template has {tmpl_n_atoms} atoms, "
                f"source array has {coords.shape[0]}."
            )

        if source_elements is not None:
            tmpl_elements = [a.GetSymbol().upper() for a in tmpl.GetAtoms()]
            src_elem_clean = [str(e).strip().upper() for e in source_elements]
            if tmpl_elements != src_elem_clean:
                raise ValueError(
                    f"FATAL: Element sequence mismatch between template and MD atoms!\n"
                    f"Template: {tmpl_elements}\nSource:   {src_elem_clean}"
                )

        # Transfer all 3D coordinates (heavy atoms and all protons) in-memory
        for i in range(tmpl_n_atoms):
            p = coords[i]
            if abs(p[0]) < 1e-4 and abs(p[1]) < 1e-4 and abs(p[2]) < 1e-4:
                raise ValueError(f"FATAL: Degenerate coordinate (0, 0, 0) detected at atom index {i}!")
            res_conf.SetAtomPosition(i, [float(p[0]), float(p[1]), float(p[2])])

    # -------------------------------------------------------------------------
    # MODALITY 2: RDKit Chem.Mol Conformation (Docking engines / GNINA / Glide)
    # -------------------------------------------------------------------------
    elif isinstance(src, Chem.Mol):
        src_mol = src
        src_conf = src_mol.GetConformer()
        src_n_atoms = src_mol.GetNumAtoms()

        tmpl_elements = [a.GetSymbol().upper() for a in tmpl.GetAtoms()]
        src_elements = [a.GetSymbol().upper() for a in src_mol.GetAtoms()]
        direct_possible = (src_n_atoms == tmpl_n_atoms and tmpl_elements == src_elements)

        if direct_possible and not relax_missing_hydrogens:
            for i in range(tmpl_n_atoms):
                pos = src_conf.GetAtomPosition(i)
                res_conf.SetAtomPosition(i, pos)
        else:
            tmpl_heavy = Chem.RemoveHs(tmpl, sanitize=False)
            src_heavy = Chem.RemoveHs(src_mol, sanitize=False)

            n_tmpl_h = tmpl_heavy.GetNumAtoms()
            n_src_h = src_heavy.GetNumAtoms()
            if n_tmpl_h != n_src_h:
                raise ValueError(f"FATAL: Heavy atom count mismatch: template={n_tmpl_h}, source={n_src_h}")

            Chem.AssignStereochemistry(tmpl_heavy, force=True, cleanIt=True)
            chiral_centers = []
            tmpl_h_conf = tmpl_heavy.GetConformer()
            for a in tmpl_heavy.GetAtoms():
                if a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
                    nbrs = [nbr.GetIdx() for nbr in a.GetNeighbors()]
                    if len(nbrs) >= 3:
                        vol = signed_chiral_volume(tmpl_h_conf, a.GetIdx(), nbrs)
                        if abs(vol) > 0.05:
                            chiral_centers.append((a.GetIdx(), nbrs, np.sign(vol)))

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
                raise ValueError("FATAL: Topological isomorphism failed between template and source molecule!")

            src_h_conf = src_heavy.GetConformer()
            best_match = None
            best_sq_dist = float("inf")

            for match in matches:
                chiral_valid = True
                for tmpl_c, tmpl_nbrs, expected_sign in chiral_centers:
                    src_c = match[tmpl_c]
                    src_nbrs = [match[n] for n in tmpl_nbrs]
                    src_vol = signed_chiral_volume(src_h_conf, src_c, src_nbrs)
                    if np.sign(src_vol) != expected_sign and abs(src_vol) > 0.05:
                        chiral_valid = False
                        break
                if not chiral_valid:
                    continue

                sq_dist = 0.0
                for t_i, s_i in enumerate(match):
                    pt = tmpl_h_conf.GetAtomPosition(t_i)
                    ps = src_h_conf.GetAtomPosition(s_i)
                    sq_dist += (pt.x - ps.x)**2 + (pt.y - ps.y)**2 + (pt.z - ps.z)**2

                if sq_dist < best_sq_dist:
                    best_sq_dist = sq_dist
                    best_match = match

            if best_match is None:
                raise ValueError("FATAL: Stereochemical parity check failed! Inversion of chiral center detected.")

            heavy_map = {}
            for t_idx, s_idx in enumerate(best_match):
                orig_t = tmpl_heavy.GetAtomWithIdx(t_idx).GetIdx()
                orig_s = src_heavy.GetAtomWithIdx(s_idx).GetIdx()
                heavy_map[orig_t] = orig_s
                pos = src_conf.GetAtomPosition(orig_s)
                res_conf.SetAtomPosition(orig_t, pos)

            if src_n_atoms < tmpl_n_atoms and relax_missing_hydrogens:
                temp_h = Chem.RemoveHs(res_mol)
                temp_h = Chem.AddHs(temp_h, addCoords=True)
                mp = AllChem.MMFFGetMoleculeProperties(temp_h)
                ff = AllChem.MMFFGetMoleculeForceField(temp_h, mp) if mp else AllChem.UFFGetMoleculeForceField(temp_h)
                if ff:
                    for a in temp_h.GetAtoms():
                        if a.GetAtomicNum() != 1:
                            ff.AddFixedPoint(a.GetIdx())
                    ff.Minimize(maxIts=500)

                temp_c = temp_h.GetConformer()
                for orig_t, orig_s in heavy_map.items():
                    tmpl_hs = [nbr.GetIdx() for nbr in tmpl.GetAtomWithIdx(orig_t).GetNeighbors() if nbr.GetAtomicNum() == 1]
                    temp_hs = [nbr.GetIdx() for nbr in temp_h.GetAtomWithIdx(orig_t).GetNeighbors() if nbr.GetAtomicNum() == 1]
                    for k in range(min(len(tmpl_hs), len(temp_hs))):
                        res_conf.SetAtomPosition(tmpl_hs[k], temp_c.GetAtomPosition(temp_hs[k]))

            for prop in src_mol.GetPropNames():
                res_mol.SetProp(prop, src_mol.GetProp(prop))
            if src_mol.HasProp("_Name"):
                res_mol.SetProp("_Name", src_mol.GetProp("_Name"))

    else:
        raise TypeError(f"FATAL: Unsupported source type: {type(src)}. Expected Chem.Mol or numpy array.")

    # -------------------------------------------------------------------------
    # FAIL-FAST PHYSICAL BARRIER: Strict Validation of Covalent Bond Lengths
    # -------------------------------------------------------------------------
    for bond in res_mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        p1 = res_conf.GetAtomPosition(a1.GetIdx())
        p2 = res_conf.GetAtomPosition(a2.GetIdx())
        dist = float(np.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2))

        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            if dist > 1.35 or dist < 0.80:
                h_a = a1 if a1.GetAtomicNum() == 1 else a2
                hv_a = a2 if a1.GetAtomicNum() == 1 else a1
                raise ValueError(
                    f"FATAL PHYSICAL VIOLATION: Severed covalent bond detected during coordinate transfer! "
                    f"Bond {hv_a.GetSymbol()}{hv_a.GetIdx()+1}-{h_a.GetSymbol()}{h_a.GetIdx()+1} "
                    f"has unphysical length {dist:.3f} Å (allowed range: 0.80 - 1.35 Å)."
                )
        else:
            if dist > 2.20 or dist < 1.00:
                raise ValueError(
                    f"FATAL PHYSICAL VIOLATION: Distorted heavy covalent bond detected! "
                    f"Bond {a1.GetSymbol()}{a1.GetIdx()+1}-{a2.GetSymbol()}{a2.GetIdx()+1} "
                    f"has unphysical length {dist:.3f} Å (allowed range: 1.00 - 2.20 Å)."
                )

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
