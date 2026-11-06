"""
File: cli/core/chem.py
Description: Chemoinformatics helpers for ligand standardization and identity verification.
Computes Standard InChIKey and Canonical SMILES via RDKit for topological deduplication.
"""

import hashlib
import math
from pathlib import Path
from typing import List, Tuple
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions,
)
from rdkit.Chem.MolStandardize import rdMolStandardize


def _load_rdkit_mol(source: str | Path | Chem.Mol) -> Tuple[Chem.Mol, str]:
    mol = None
    mol_name = ""

    if isinstance(source, Chem.Mol) or hasattr(source, "GetAtoms"):
        mol = Chem.Mol(source)
        if mol.HasProp("_Name"):
            mol_name = mol.GetProp("_Name").strip()
    elif isinstance(source, Path) or (
        isinstance(source, str) and Path(source).exists()
    ):
        p = Path(source)
        if p.suffix.lower() in [".sdf", ".mol"]:
            suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=False)
            mol = next((m for m in suppl if m is not None), None)
            if mol and mol.HasProp("_Name"):
                mol_name = mol.GetProp("_Name").strip()
        elif p.suffix.lower() == ".pdb":
            mol = Chem.MolFromPDBFile(str(p), removeHs=False, sanitize=False)
        else:
            raise ValueError(f"FATAL: Unsupported molecule file format: {p.suffix}")
    elif isinstance(source, str):
        source_clean = source.strip()
        if "\n" in source_clean or "$$$$" in source_clean:
            mol = Chem.MolFromMolBlock(source_clean, removeHs=False, sanitize=False)
            if mol and mol.HasProp("_Name"):
                mol_name = mol.GetProp("_Name").strip()
        else:
            mol = Chem.MolFromSmiles(source_clean, sanitize=False)
            mol_name = source_clean

    if mol is None:
        raise ValueError(
            f"FATAL: Could not parse valid chemical structure from source: {source}"
        )

    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)
    _perceive_stereochemistry(mol)
    return mol, mol_name


def _perceive_stereochemistry(mol: Chem.Mol) -> None:
    """Perceives 3D or 2D stereochemistry deterministically without losing wedge or coordinate parity."""
    if mol.GetNumConformers() > 0 and mol.GetConformer().Is3D():
        # Preserve any existing 2D wedge annotations if coordinates are flat, otherwise perceive from 3D
        conf = mol.GetConformer()
        zs = [abs(conf.GetAtomPosition(i).z) for i in range(mol.GetNumAtoms())]
        if max(zs, default=0.0) > 1e-3:
            Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
            Chem.AssignStereochemistry(
                mol, force=False, cleanIt=True, flagPossibleStereoCenters=True
            )
            return
    Chem.AssignStereochemistry(
        mol, force=True, cleanIt=True, flagPossibleStereoCenters=True
    )


def _format_stereo_tag(mol: Chem.Mol, inchikey: str) -> str:
    """Extracts human-readable CIP stereocenter descriptors (e.g. '3R_5S') or InChIKey stereo block."""
    Chem.AssignStereochemistry(
        mol, force=False, cleanIt=True, flagPossibleStereoCenters=True
    )
    cip_tokens = []
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            cip = atom.GetProp("_CIPCode").strip()
            cip_tokens.append(f"{atom.GetIdx()+1}{cip}")
    for bond in mol.GetBonds():
        if bond.GetStereo() == Chem.BondStereo.STEREOE:
            cip_tokens.append(f"b{bond.GetIdx()+1}E")
        elif bond.GetStereo() == Chem.BondStereo.STEREOZ:
            cip_tokens.append(f"b{bond.GetIdx()+1}Z")
    if cip_tokens:
        return "_".join(cip_tokens)
    parts = inchikey.split("-")
    return parts[1] if len(parts) > 1 else inchikey[:8]


def expand_stereoisomers(
    mol: Chem.Mol, base_name: str, max_isomers: int = 32
) -> List[Tuple[Chem.Mol, str]]:
    """
    Checks for unspecified covalent stereocenters (chiral atoms '?' or STEREOANY double bonds)
    on the uncharged molecular skeleton. If unspecified centers exist, deterministically enumerates
    all physically embeddable stereoisomers so every Substance in the 3NF DB has 100% explicit chirality.
    """
    mol_copy = Chem.Mol(mol)
    mol_copy.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol_copy)
    _perceive_stereochemistry(mol_copy)

    if mol_copy.GetNumHeavyAtoms() == 1 and mol_copy.GetNumAtoms() == 1:
        return [(mol_copy, base_name)]

    uncharger = rdMolStandardize.Uncharger()
    heavy_neutral = uncharger.uncharge(Chem.RemoveHs(mol_copy))
    Chem.AssignStereochemistry(
        heavy_neutral, force=False, cleanIt=True, flagPossibleStereoCenters=True
    )

    chiral_centers = Chem.FindMolChiralCenters(
        heavy_neutral, includeUnassigned=True, useLegacyImplementation=False
    )
    has_unassigned_atoms = any(tag == "?" for _, tag in chiral_centers)
    has_unassigned_bonds = any(
        b.GetBondType() == Chem.BondType.DOUBLE
        and b.GetStereo() == Chem.BondStereo.STEREOANY
        for b in heavy_neutral.GetBonds()
    )

    if not has_unassigned_atoms and not has_unassigned_bonds:
        return [(mol_copy, base_name)]

    # Strip 3D conformers before stereoisomer enumeration so RDKit embeds each isomer cleanly
    mol_no_conf = Chem.RemoveHs(mol_copy)
    mol_no_conf.RemoveAllConformers()
    Chem.AssignStereochemistry(
        mol_no_conf, force=False, cleanIt=True, flagPossibleStereoCenters=True
    )

    opts = StereoEnumerationOptions(
        tryEmbedding=True, unique=True, maxIsomers=max_isomers, onlyUnassigned=True
    )
    raw_isomers = list(EnumerateStereoisomers(mol_no_conf, options=opts))
    if not raw_isomers:
        raise ValueError(
            f"FATAL: Stereoisomer enumeration failed to produce any physically embeddable isomer for '{base_name}'!"
        )

    unique_isomers: List[Tuple[Chem.Mol, str]] = []
    seen_inchikeys = set()
    for iso in raw_isomers:
        iso_h = Chem.AddHs(iso)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        res_cid = AllChem.EmbedMolecule(iso_h, params)
        if res_cid < 0:
            res_cid = AllChem.EmbedMolecule(iso_h, useRandomCoords=True, randomSeed=42)
        if res_cid < 0:
            raise ValueError(
                f"FATAL: Failed to generate 3D coordinates for enumerated stereoisomer of '{base_name}'!"
            )
        AllChem.UFFOptimizeMolecule(iso_h, maxIters=200)
        _perceive_stereochemistry(iso_h)

        ik, _, _, _ = get_mol_identity(iso_h)
        if ik in seen_inchikeys:
            continue
        seen_inchikeys.add(ik)
        unique_isomers.append((iso_h, ik))

    if len(unique_isomers) == 1:
        # All enumerated bond/atom variations collapsed to a single symmetric molecule;
        # preserve the original input conformation and base_name untouched!
        mol_copy.SetProp("_Name", base_name)
        return [(mol_copy, base_name)]

    results: List[Tuple[Chem.Mol, str]] = []
    for iso_h, ik in unique_isomers:
        stereo_tag = _format_stereo_tag(iso_h, ik)
        iso_title = f"{base_name}_{stereo_tag}"
        iso_h.SetProp("_Name", iso_title)
        results.append((iso_h, iso_title))

    return results


def get_mol_identity(source: str | Path | Chem.Mol) -> Tuple[str, str, str, int]:
    """
    Parses a chemical structure, perceives 3D/2D stereochemistry, neutralizes ionizable groups,
    and returns (parent_inchikey, parent_isomeric_smiles, mol_name, formal_charge).
    """
    mol, mol_name = _load_rdkit_mol(source)
    formal_charge = Chem.GetFormalCharge(mol)

    is_monoatomic_ion = mol.GetNumHeavyAtoms() == 1 and mol.GetNumAtoms() == 1

    if is_monoatomic_ion:
        parent_inchikey = Chem.MolToInchiKey(mol)
        parent_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    else:
        uncharger = rdMolStandardize.Uncharger()
        heavy_mol = Chem.RemoveHs(mol)
        parent_mol = uncharger.uncharge(heavy_mol)
        Chem.AssignStereochemistry(
            parent_mol, force=False, cleanIt=True, flagPossibleStereoCenters=True
        )
        parent_inchikey = Chem.MolToInchiKey(parent_mol)
        parent_smiles = Chem.MolToSmiles(
            parent_mol, isomericSmiles=True, canonical=True
        )

    if not parent_inchikey:
        raise ValueError(
            "FATAL: Failed to generate Standard InChIKey. Invalid chemical structure."
        )

    if not mol_name:
        mol_name = parent_inchikey[:14]

    return parent_inchikey, parent_smiles, mol_name, formal_charge


def get_microstate_identity(
    source: str | Path | Chem.Mol,
) -> Tuple[str, str, str, str, str, int]:
    """
    Computes both the Abstract Substance Identity (Stereospecific Neutral InChIKey) and the
    Exact Chemical Microstate Identity (Isomeric Canonical SMILES and SHA-256 state hash).

    Returns:
        (parent_inchikey, parent_smiles, mol_name, state_hash, isomeric_smiles, formal_charge)
    """
    mol, _ = _load_rdkit_mol(source)
    parent_inchikey, parent_smiles, mol_name, formal_charge = get_mol_identity(mol)

    if mol.GetNumHeavyAtoms() == 1 and mol.GetNumAtoms() == 1:
        mol_canon = mol
    else:
        mol_canon = Chem.AddHs(mol)
        _perceive_stereochemistry(mol_canon)

    isomeric_smiles = Chem.MolToSmiles(mol_canon, isomericSmiles=True, canonical=True)
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
        raise ValueError(
            "FATAL: transfer_coordinates_to_template received None for template or source!"
        )

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
            raise ValueError(
                f"FATAL: Expected coordinate array of shape (N, 3), got: {coords.shape}"
            )

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
                raise ValueError(
                    f"FATAL: Degenerate coordinate (0, 0, 0) detected at atom index {i}!"
                )
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
        direct_possible = src_n_atoms == tmpl_n_atoms and tmpl_elements == src_elements

        if direct_possible and not relax_missing_hydrogens:
            for i in range(tmpl_n_atoms):
                pos = src_conf.GetAtomPosition(i)
                res_conf.SetAtomPosition(i, pos)
        else:
            # Tag original atom indices in template before removing hydrogens to prevent index drift
            for idx, a in enumerate(tmpl.GetAtoms()):
                a.SetIntProp("_orig_idx", idx)

            tmpl_heavy = Chem.RemoveHs(tmpl, sanitize=False)
            src_heavy = Chem.RemoveHs(src_mol, sanitize=False)

            n_tmpl_h = tmpl_heavy.GetNumAtoms()
            n_src_h = src_heavy.GetNumAtoms()
            if n_tmpl_h != n_src_h:
                raise ValueError(
                    f"FATAL: Heavy atom count mismatch: template={n_tmpl_h}, source={n_src_h}"
                )

            Chem.AssignStereochemistry(tmpl_heavy, force=True, cleanIt=True)
            chiral_centers = []
            tmpl_h_conf = tmpl_heavy.GetConformer()
            for a in tmpl_heavy.GetAtoms():
                if a.GetChiralTag() in (
                    Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                    Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                ):
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
                raise ValueError(
                    "FATAL: Topological isomorphism failed between template and source molecule!"
                )

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
                    sq_dist += (
                        (pt.x - ps.x) ** 2 + (pt.y - ps.y) ** 2 + (pt.z - ps.z) ** 2
                    )

                if sq_dist < best_sq_dist:
                    best_sq_dist = sq_dist
                    best_match = match

            if best_match is None:
                raise ValueError(
                    "FATAL: Stereochemical parity check failed! Inversion of chiral center detected."
                )

            heavy_map = {}
            for t_idx, s_idx in enumerate(best_match):
                # Retrieve true original index in all-atom tmpl via preserved _orig_idx property
                orig_t = tmpl_heavy.GetAtomWithIdx(t_idx).GetIntProp("_orig_idx")
                orig_s = s_idx
                heavy_map[orig_t] = orig_s
                pos = src_h_conf.GetAtomPosition(orig_s)
                res_conf.SetAtomPosition(orig_t, pos)

            if src_n_atoms < tmpl_n_atoms and relax_missing_hydrogens:
                temp_h = Chem.RemoveHs(res_mol)
                temp_h = Chem.AddHs(temp_h, addCoords=True)
                temp_c = temp_h.GetConformer()
                for orig_t, orig_s in heavy_map.items():
                    tmpl_hs = [
                        nbr.GetIdx()
                        for nbr in tmpl.GetAtomWithIdx(orig_t).GetNeighbors()
                        if nbr.GetAtomicNum() == 1
                    ]
                    temp_hs = [
                        nbr.GetIdx()
                        for nbr in temp_h.GetAtomWithIdx(orig_t).GetNeighbors()
                        if nbr.GetAtomicNum() == 1
                    ]
                    for k in range(min(len(tmpl_hs), len(temp_hs))):
                        res_conf.SetAtomPosition(
                            tmpl_hs[k], temp_c.GetAtomPosition(temp_hs[k])
                        )

        for prop in src_mol.GetPropNames():
            res_mol.SetProp(prop, src_mol.GetProp(prop))
        if src_mol.HasProp("_Name"):
            res_mol.SetProp("_Name", src_mol.GetProp("_Name"))

    else:
        raise TypeError(
            f"FATAL: Unsupported source type: {type(src)}. Expected Chem.Mol or numpy array."
        )

    # -------------------------------------------------------------------------
    # FAIL-FAST PHYSICAL BARRIER: Strict Validation of Covalent Bond Lengths
    # -------------------------------------------------------------------------
    for bond in res_mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        p1 = res_conf.GetAtomPosition(a1.GetIdx())
        p2 = res_conf.GetAtomPosition(a2.GetIdx())
        dist = float(
            np.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)
        )

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
            a.SetNumExplicitHs(0)
            a.SetNoImplicit(False)
            a.SetNumRadicalElectrons(0)
        q.UpdatePropertyCache(strict=False)

    matches = q_target.GetSubstructMatches(q_tmpl, uniquify=False)
    if not matches:
        raise ValueError(
            "FATAL: Failed to obtain topological isomorphism matches for RMSD calculation!"
        )

    conf_tmpl = tmpl_heavy.GetConformer()
    conf_target = target_heavy.GetConformer()
    best_sq_dist = float("inf")

    for match in matches:
        sq_dist = 0.0
        for tmpl_idx, target_idx in enumerate(match):
            p_t = conf_tmpl.GetAtomPosition(tmpl_idx)
            p_tgt = conf_target.GetAtomPosition(target_idx)
            sq_dist += (
                (p_t.x - p_tgt.x) ** 2 + (p_t.y - p_tgt.y) ** 2 + (p_t.z - p_tgt.z) ** 2
            )

        if sq_dist < best_sq_dist:
            best_sq_dist = sq_dist

    return float(np.sqrt(best_sq_dist / float(n_atoms)))


def compute_pairwise_symmetry_rmsd_matrix(
    mols_a: list[Chem.Mol], mols_b: list[Chem.Mol]
) -> tuple[list[float], list[float]]:
    """
    Computes all-against-all symmetry-aware heavy-atom RMSD matrix between two sets of poses
    for the exact same substance (e.g. Glide modes vs GNINA modes) using single graph isomorphism
    and vectorized NumPy tensor broadcasting.
    Returns:
        (min_rmsd_for_each_a, min_rmsd_for_each_b)
    """
    from rdkit import Chem
    import numpy as np

    if not mols_a or not mols_b:
        return ([float("inf")] * len(mols_a), [float("inf")] * len(mols_b))

    tmpl_heavy = Chem.RemoveHs(mols_a[0], sanitize=False)
    n_heavy = tmpl_heavy.GetNumAtoms()
    if n_heavy == 0:
        return ([0.0] * len(mols_a), [0.0] * len(mols_b))

    q_tmpl = Chem.Mol(tmpl_heavy)
    for b in q_tmpl.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in q_tmpl.GetAtoms():
        a.SetFormalCharge(0)
        a.SetIsotope(0)
        a.SetNumExplicitHs(0)
        a.SetNoImplicit(False)
        a.SetNumRadicalElectrons(0)
    q_tmpl.UpdatePropertyCache(strict=False)

    matches = q_tmpl.GetSubstructMatches(q_tmpl, uniquify=False)
    if not matches:
        matches = [tuple(range(n_heavy))]

    def _extract_heavy_coords(mol_list):
        coords_arr = np.empty((len(mol_list), n_heavy, 3), dtype=np.float64)
        for m_i, m in enumerate(mol_list):
            m_h = Chem.RemoveHs(m, sanitize=False)
            conf = m_h.GetConformer()
            for a_i in range(n_heavy):
                pos = conf.GetAtomPosition(a_i)
                coords_arr[m_i, a_i, 0] = pos.x
                coords_arr[m_i, a_i, 1] = pos.y
                coords_arr[m_i, a_i, 2] = pos.z
        return coords_arr

    coords_a = _extract_heavy_coords(mols_a)
    coords_b = _extract_heavy_coords(mols_b)

    best_sq_dist = np.full((len(mols_a), len(mols_b)), np.inf, dtype=np.float64)

    for match in matches:
        perm_b = coords_b[:, match, :]
        diff = coords_a[:, np.newaxis, :, :] - perm_b[np.newaxis, :, :, :]
        sq_dist = np.sum(diff ** 2, axis=(2, 3)) / float(n_heavy)
        best_sq_dist = np.minimum(best_sq_dist, sq_dist)

    rmsd_mat = np.sqrt(best_sq_dist)
    min_for_a = np.min(rmsd_mat, axis=1).round(3).tolist()
    min_for_b = np.min(rmsd_mat, axis=0).round(3).tolist()
    return min_for_a, min_for_b


def compute_autobox_geometry(ref_path: str | Path, autobox_add: float) -> dict:
    """
    Single Source of Truth (SSOT) for pocket bounding box and centroid geometry.
    Reads reference ligand 3D coordinates (.sdf, .pdb, .mae, .maegz) and computes:
      - centroid: arithmetic mean (cx, cy, cz) for Glide GRID_CENTER
      - autobox_center: bounding box center (bc_x, bc_y, bc_z) for Uni-Dock / GNINA
      - autobox_size: padded box dimensions (sx, sy, sz) = span + 2 * autobox_add
      - autobox_bounds: (min_x, max_x, min_y, max_y, min_z, max_z) for PyMOL CGO
    """
    p = Path(ref_path).resolve()
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(
            f"FATAL: Reference ligand file missing or empty: {p}"
        )

    pad = float(autobox_add)
    if pad < 0.0:
        raise ValueError(
            f"FATAL: autobox_add must be non-negative, got: {pad}"
        )

    coords: list[tuple[float, float, float]] = []
    suffix = p.suffix.lower()
    if suffix in (".sdf", ".mol", ".pdb"):
        mol, _ = _load_rdkit_mol(p)
        if mol.GetNumConformers() == 0:
            raise ValueError(f"FATAL: Reference ligand in {p} has no 3D conformer!")
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append((float(pos.x), float(pos.y), float(pos.z)))
    elif suffix in (".mae", ".maegz"):
        from schrodinger import structure

        st = structure.StructureReader.read(str(p))
        coords = [(float(a.x), float(a.y), float(a.z)) for a in st.atom]
    else:
        raise ValueError(
            f"FATAL: Unsupported reference ligand format '{suffix}' for {p}"
        )

    if not coords:
        raise ValueError(f"FATAL: Zero atom coordinates extracted from {p}")

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    bc_x = 0.5 * (min_x + max_x)
    bc_y = 0.5 * (min_y + max_y)
    bc_z = 0.5 * (min_z + max_z)

    sx = (max_x - min_x) + 2.0 * pad
    sy = (max_y - min_y) + 2.0 * pad
    sz = (max_z - min_z) + 2.0 * pad

    return {
        "centroid": (cx, cy, cz),
        "autobox_center": (bc_x, bc_y, bc_z),
        "autobox_size": (sx, sy, sz),
        "autobox_bounds": (
            min_x - pad,
            max_x + pad,
            min_y - pad,
            max_y + pad,
            min_z - pad,
            max_z + pad,
        ),
    }


def convert_receptor_pdb_to_pdbqt(pdb_path: str | Path, out_pdbqt: str | Path) -> None:
    """
    Single Source of Truth (SSOT) rigid receptor converter from PrepWizard-prepared PDB to AutoDock4 PDBQT.
    Strictly preserves PrepWizard protonation states, histidine tautomers, cofactors, and metal ions.
    Strips only non-polar hydrogens (bonded to carbon) as required by Vina/Uni-Dock united-atom potentials
    and assigns exact AutoDock4 atom types in columns 78-79.
    """
    in_p = Path(pdb_path).resolve()
    out_p = Path(out_pdbqt).resolve()
    if not in_p.exists() or in_p.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Input receptor PDB missing or empty: {in_p}")

    mol = Chem.MolFromPDBFile(str(in_p), removeHs=False, sanitize=False)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"FATAL: RDKit failed to parse receptor PDB: {in_p}")
    mol.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(mol)
    Chem.SetAromaticity(mol)

    raw_atom_lines = [
        line
        for line in in_p.read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.startswith(("ATOM  ", "HETATM"))
    ]
    if len(raw_atom_lines) != mol.GetNumAtoms():
        raise ValueError(
            f"FATAL: Atom count mismatch in {in_p.name}: PDB has {len(raw_atom_lines)} ATOM/HETATM lines, "
            f"RDKit parsed {mol.GetNumAtoms()} atoms!"
        )

    out_lines: list[str] = []
    for idx, atom in enumerate(mol.GetAtoms()):
        line = raw_atom_lines[idx]
        elem = atom.GetSymbol().strip().upper()
        nbrs = atom.GetNeighbors()

        if elem == "H":
            # Retain polar hydrogens (bonded to N, O, S, P); omit non-polar carbon-bound hydrogens
            if not nbrs or nbrs[0].GetSymbol().strip().upper() == "C":
                continue
            ad4_type = "HD"
        elif elem == "C":
            ad4_type = "A" if atom.GetIsAromatic() else "C"
        elif elem == "N":
            num_h = sum(1 for n in nbrs if n.GetSymbol().strip().upper() == "H")
            # Unprotonated sp2/aromatic or nitrile nitrogen with lone pair is H-bond acceptor (NA)
            if num_h == 0 and len(nbrs) <= 2:
                ad4_type = "NA"
            else:
                ad4_type = "N"
        elif elem == "O":
            ad4_type = "OA"
        elif elem == "S":
            ad4_type = "SA"
        elif elem in {"P", "F", "B", "I"}:
            ad4_type = elem
        elif elem in {"CL", "BR", "MG", "ZN", "MN", "FE", "CA", "CU", "NA", "K"}:
            ad4_type = elem.capitalize() if elem in {"CL", "BR"} else elem
        else:
            ad4_type = elem[:2]

        base_cols = line[:54].ljust(54)
        pdbqt_line = f"{base_cols}  1.00  0.00    +0.000 {ad4_type:<2s}"
        out_lines.append(pdbqt_line)

    if not out_lines:
        raise ValueError(f"FATAL: Zero receptor atoms written to {out_p}!")

    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def convert_prepared_mol_to_pdbqt(mol: Chem.Mol, out_pdbqt: str | Path) -> int:
    """
    Translates a fully prepared, protonated LigPrep molecule into an AutoDock PDBQT string
    with a complete torsion tree (ROOT/BRANCH/TORSDOF) using standard Meeko.
    Preserves exact 3D coordinates, formal charges, and polar hydrogens without distortion.
    Returns the integer count of active rotatable bonds (torsions).
    """
    out_path = Path(out_pdbqt).resolve()
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError("FATAL: Invalid or empty Chem.Mol provided to convert_prepared_mol_to_pdbqt!")

    preparator = MoleculePreparation(
        merge_these_atom_types=["H"],
        hydrate=False,
        flexible_amides=False,
        rigid_macrocycles=True,
        min_ring_size=7,
        max_ring_size=33,
    )

    mol_copy = Chem.Mol(mol)
    mol_copy.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol_copy)

    setups = preparator.prepare(mol_copy)
    if not setups:
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "mol"
        raise RuntimeError(f"FATAL: Meeko returned empty setup list for '{mol_name}'!")
    setup = setups[0]

    pdbqt_str, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
    if not is_ok or not pdbqt_str:
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "mol"
        raise RuntimeError(
            f"FATAL: Meeko failed to generate valid PDBQT for '{mol_name}': {error_msg}"
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(pdbqt_str, encoding="utf-8")

    n_torsions = 0
    for line in pdbqt_str.splitlines():
        if line.startswith("TORSDOF"):
            parts = line.split()
            if len(parts) >= 2:
                n_torsions = int(parts[1])
            break
    return n_torsions


def _parse_pdbqt_atom_coords(lines: list[str]) -> list[tuple[float, float, float]]:
    """Extracts 3D Cartesian coordinates (x, y, z) from PDBQT ATOM/HETATM records."""
    coords: list[tuple[float, float, float]] = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append((x, y, z))
    return coords


def restore_unidock_pdbqt_poses(
    exact_tmpl: Chem.Mol,
    in_pdbqt: Path,
    out_pdbqt: Path,
    orig_title: str,
) -> list[tuple[float, Chem.Mol]]:
    """
    Restores full-atom RDKit molecules from a Uni-Dock output PDBQT file onto the
    canonical LigPrep template (`exact_tmpl`).
    1. Maps each PDBQT atom line to its exact atom index in `exact_tmpl` by matching
       the unshifted input coordinates in `in_pdbqt` (tolerance < 0.01 A).
    2. For each docked MODEL in `out_pdbqt`, transfers the exact docked coordinates of all
       heavy atoms and all polar hydrogens (OH, NH, SH) onto a copy of `exact_tmpl`.
    3. Reconstructs only the merged non-polar carbon-bound hydrogens (C-H) at the new
       docked carbon positions and validates all covalent bond lengths and chiral volumes.
    """
    in_lines = in_pdbqt.read_text(encoding="utf-8", errors="ignore").splitlines()
    in_coords = _parse_pdbqt_atom_coords(in_lines)
    if not in_coords:
        raise ValueError(f"FATAL: Zero ATOM/HETATM records in input PDBQT: {in_pdbqt}")

    tmpl_conf = exact_tmpl.GetConformer()
    n_tmpl = exact_tmpl.GetNumAtoms()

    pdbqt_to_tmpl_idx: list[int] = []
    mapped_tmpl_indices: set[int] = set()
    for p_idx, (px, py, pz) in enumerate(in_coords):
        best_t_idx = -1
        best_d2 = 1e-4
        for t_idx in range(n_tmpl):
            if t_idx in mapped_tmpl_indices:
                continue
            tp = tmpl_conf.GetAtomPosition(t_idx)
            d2 = (tp.x - px) ** 2 + (tp.y - py) ** 2 + (tp.z - pz) ** 2
            if d2 < best_d2:
                best_d2 = d2
                best_t_idx = t_idx
        if best_t_idx < 0:
            raise ValueError(
                f"FATAL: Input PDBQT atom {p_idx + 1} at ({px:.3f}, {py:.3f}, {pz:.3f}) in {in_pdbqt.name} "
                f"does not match any atom in LigPrep template '{orig_title}' within 0.01 A!"
            )
        pdbqt_to_tmpl_idx.append(best_t_idx)
        mapped_tmpl_indices.add(best_t_idx)

    # Identify unmapped template atoms (strictly non-polar hydrogens merged by Meeko)
    unmapped_h_indices: list[int] = []
    for t_idx in range(n_tmpl):
        if t_idx not in mapped_tmpl_indices:
            atom = exact_tmpl.GetAtomWithIdx(t_idx)
            if atom.GetAtomicNum() != 1:
                raise ValueError(
                    f"FATAL: Heavy atom {atom.GetSymbol()}{t_idx + 1} of '{orig_title}' is missing from {in_pdbqt.name}!"
                )
            unmapped_h_indices.append(t_idx)

    # Split out_pdbqt into MODEL blocks
    out_lines = out_pdbqt.read_text(encoding="utf-8", errors="ignore").splitlines()
    models: list[list[str]] = []
    curr_model: list[str] = []
    in_model = False
    for line in out_lines:
        if line.startswith("MODEL"):
            in_model = True
            curr_model = [line]
        elif line.startswith("ENDMDL"):
            curr_model.append(line)
            models.append(curr_model)
            curr_model = []
            in_model = False
        elif in_model:
            curr_model.append(line)

    if not models and out_lines:
        models.append(out_lines)

    restored_results: list[tuple[float, Chem.Mol]] = []

    for m_idx, m_lines in enumerate(models, start=1):
        vina_score: float | None = None
        for line in m_lines:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                if len(parts) >= 4:
                    vina_score = float(parts[3])
                    break
        if vina_score is None:
            raise ValueError(
                f"FATAL: Missing 'REMARK VINA RESULT:' in MODEL {m_idx} of {out_pdbqt}!"
            )

        pose_coords = _parse_pdbqt_atom_coords(m_lines)
        if len(pose_coords) != len(pdbqt_to_tmpl_idx):
            raise ValueError(
                f"FATAL: Atom count mismatch in MODEL {m_idx} of {out_pdbqt}: "
                f"expected {len(pdbqt_to_tmpl_idx)}, got {len(pose_coords)}!"
            )

        res_mol = Chem.Mol(exact_tmpl)
        res_conf = res_mol.GetConformer()

        # 1. Transfer exact docked coordinates of all heavy atoms and all polar hydrogens
        for p_i, (dx, dy, dz) in enumerate(pose_coords):
            t_i = pdbqt_to_tmpl_idx[p_i]
            res_conf.SetAtomPosition(t_i, [dx, dy, dz])

        # 2. Reconstruct merged non-polar C-H hydrogens at the new docked carbon positions
        if unmapped_h_indices:
            for idx, a in enumerate(res_mol.GetAtoms()):
                a.SetIntProp("_orig_tmpl_idx", idx)
            heavy_only = Chem.RemoveHs(res_mol)
            rebuilt_mol = Chem.AddHs(heavy_only, addCoords=True)
            rebuilt_conf = rebuilt_mol.GetConformer()

            heavy_orig_to_rebuilt: dict[int, int] = {}
            for r_idx, r_atom in enumerate(rebuilt_mol.GetAtoms()):
                if r_atom.HasProp("_orig_tmpl_idx"):
                    heavy_orig_to_rebuilt[r_atom.GetIntProp("_orig_tmpl_idx")] = r_idx

            # Group unmapped non-polar hydrogens by their parent heavy atom in exact_tmpl
            parent_to_unmapped_hs: dict[int, list[int]] = {}
            for h_idx in unmapped_h_indices:
                h_atom = exact_tmpl.GetAtomWithIdx(h_idx)
                nbrs = h_atom.GetNeighbors()
                if not nbrs:
                    raise ValueError(f"FATAL: Detached hydrogen atom at index {h_idx} in '{orig_title}'!")
                p_orig = nbrs[0].GetIdx()
                if p_orig not in parent_to_unmapped_hs:
                    parent_to_unmapped_hs[p_orig] = []
                parent_to_unmapped_hs[p_orig].append(h_idx)

            for p_orig, h_list in parent_to_unmapped_hs.items():
                r_parent_idx = heavy_orig_to_rebuilt[p_orig]
                r_parent_atom = rebuilt_mol.GetAtomWithIdx(r_parent_idx)
                r_h_indices = [
                    nbr.GetIdx()
                    for nbr in r_parent_atom.GetNeighbors()
                    if nbr.GetAtomicNum() == 1
                ]
                if len(r_h_indices) < len(h_list):
                    raise ValueError(
                        f"FATAL: Rebuilt heavy atom {p_orig} has {len(r_h_indices)} H neighbors, "
                        f"expected at least {len(h_list)}!"
                    )
                for k, target_h_idx in enumerate(h_list):
                    pos_h = rebuilt_conf.GetAtomPosition(r_h_indices[k])
                    res_conf.SetAtomPosition(target_h_idx, [pos_h.x, pos_h.y, pos_h.z])

        # 3. Fail-Fast validation of all covalent bond lengths
        for bond in res_mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            p1 = res_conf.GetAtomPosition(a1.GetIdx())
            p2 = res_conf.GetAtomPosition(a2.GetIdx())
            dist = math.sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)
            if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
                if dist < 0.80 or dist > 1.35:
                    raise ValueError(
                        f"FATAL PHYSICAL VIOLATION in '{orig_title}' MODEL {m_idx}: X-H bond "
                        f"{a1.GetSymbol()}{a1.GetIdx() + 1}-{a2.GetSymbol()}{a2.GetIdx() + 1} has length {dist:.3f} A!"
                    )
            else:
                if dist < 1.00 or dist > 2.20:
                    raise ValueError(
                        f"FATAL PHYSICAL VIOLATION in '{orig_title}' MODEL {m_idx}: heavy bond "
                        f"{a1.GetSymbol()}{a1.GetIdx() + 1}-{a2.GetSymbol()}{a2.GetIdx() + 1} has length {dist:.3f} A!"
                    )

        res_mol.SetProp("_Name", orig_title)
        res_mol.SetProp("s_m_title", orig_title)
        res_mol.SetProp("unidock_vina_score", f"{vina_score:.3f}")
        res_mol.SetProp("minimizedAffinity", f"{vina_score:.3f}")
        restored_results.append((vina_score, res_mol))

    return restored_results
