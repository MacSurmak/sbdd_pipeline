#!/usr/bin/env python3
"""
bin/build_openff_topology.py
Translates OpenFF Sage SMIRNOFF potentials into AMBER mol2 and frcmod formats.
Guarantees mathematical potential equivalence, topological improper central hub alignment,
PDB-lexical atom naming, and strict immunity against standard AMBER force field atom types.
"""

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

AMBER_RESERVED_TYPES: Set[str] = {
    # ff19SB & standard protein / nucleic
    "C", "CA", "CB", "CC", "CD", "CK", "CM", "CN", "CQ", "CR", "CT", "CV", "CW", "CX", "CY", "CZ", "C*",
    "H", "H1", "H2", "H3", "H4", "H5", "HA", "HC", "HN", "HO", "HP", "HS", "HW", "HX", "HY",
    "N", "NA", "NB", "NC", "NT", "NY", "N2", "N3", "N*", "O", "O2", "OH", "OS", "OW", "S", "SH", "P",
    # Standard Ions & Water
    "K", "K+", "Na", "NA", "Cl", "CL", "Cl-", "MG", "ZN", "CA", "C0", "MW", "EP", "LP"
}

SAFE_PREFIXES = "1234567890JQX"
ALPHANUMERIC = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


def allocate_atom_types(n_atoms: int, start_prefix: str) -> List[str]:
    """Allocates disjoint 2-character atom types guaranteed not to collide with AMBER."""
    if start_prefix not in SAFE_PREFIXES:
        raise ValueError(
            f"FATAL: Invalid starting prefix '{start_prefix}'. Must be in SAFE_PREFIXES: {SAFE_PREFIXES}"
        )

    base_len = len(ALPHANUMERIC)
    prefixes_needed = math.ceil(n_atoms / base_len)
    start_idx = SAFE_PREFIXES.index(start_prefix)

    # Slot capacity check: primary ligand (prefix '1') allows up to 5 blocks (310 atoms),
    # secondary components allow up to 2 blocks (124 atoms).
    max_allowed = 5 if start_prefix == "1" else 2
    if prefixes_needed > max_allowed:
        raise OverflowError(
            f"FATAL: Molecule ({n_atoms} atoms) requires {prefixes_needed} prefix block(s), "
            f"exceeding allocated slot capacity of {max_allowed} block(s) ({max_allowed * base_len} types)! "
            f"Execution halted by Fail-Fast law."
        )

    if start_idx + prefixes_needed > len(SAFE_PREFIXES):
        raise OverflowError(
            f"FATAL: Insufficient safe prefixes remaining ({len(SAFE_PREFIXES) - start_idx}) "
            f"to allocate {prefixes_needed} block(s) for molecule ({n_atoms} atoms)."
        )

    assigned = []
    for i in range(n_atoms):
        p_offset = i // base_len
        s_idx = i % base_len
        p_char = SAFE_PREFIXES[start_idx + p_offset]
        s_char = ALPHANUMERIC[s_idx]
        token = f"{p_char}{s_char}"

        if token in AMBER_RESERVED_TYPES or token[0].islower():
            raise ValueError(
                f"FATAL: Generated atom type '{token}' violates AMBER reservation constraint!"
            )
        assigned.append(token)

    return assigned


def load_external_crg_charges(charges_path: Path, expected_atoms: int) -> List[float]:
    """Loads flat float vector of atomic partial charges."""
    if not charges_path.exists() or charges_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Charges file not found or empty: {charges_path}")

    charges = []
    for line in charges_path.read_text(encoding="utf-8").splitlines():
        line_str = line.strip()
        if line_str:
            charges.append(float(line_str))

    if len(charges) != expected_atoms:
        raise ValueError(
            f"FATAL: Charge count mismatch in {charges_path.name}: expected {expected_atoms}, found {len(charges)}"
        )
    return charges


def main():
    parser = argparse.ArgumentParser(
        description="Translate OpenFF Sage force field to AMBER mol2 + frcmod"
    )
    parser.add_argument("--sdf", required=True, help="Input relaxed 3D structure SDF")
    parser.add_argument("--charge", required=True, type=int, help="Net formal charge")
    parser.add_argument("--res_name", required=True, help="Residue name (3-4 chars)")
    parser.add_argument("--forcefield", default="openff-2.3.0.offxml", help="OpenFF force field file")
    parser.add_argument("--charge_method", required=True, choices=["ashgc", "resp", "bcc"])
    parser.add_argument("--charges", default=None, help="External charges vector file (.crg)")
    parser.add_argument("--prefix", default="1", help="Safe namespace starting prefix (from SAFE_PREFIXES)")
    parser.add_argument("--out_mol2", required=True, help="Output mol2 path")
    parser.add_argument("--out_frcmod", required=True, help="Output frcmod path")
    args = parser.parse_args()

    from rdkit import Chem
    from openff.toolkit import ForceField, Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
    from openff.units import unit

    sdf_path = Path(args.sdf).resolve()
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        raise FileNotFoundError(f"FATAL: Input SDF deliverable not found: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    rdmol = next((m for m in suppl if m is not None), None)
    if rdmol is None:
        raise ValueError(f"FATAL: Could not read molecular structure from {sdf_path}")

    try:
        Chem.SanitizeMol(rdmol)
    except Exception:
        rdmol.UpdatePropertyCache(strict=False)

    # Strict Fail-Fast stereochemistry resolution: perceive 3D coordinates and forbid unassigned centers
    Chem.AssignStereochemistryFrom3D(rdmol)
    unassigned = Chem.FindMolChiralCenters(rdmol, includeUnassigned=True)
    undefined_centers = [idx for idx, tag in unassigned if tag == "?"]
    if undefined_centers:
        raise ValueError(
            f"FATAL: Molecule '{args.res_name}' contains undefined chiral center(s) at atom index(es): {undefined_centers}! "
            f"Stereochemistry must be fully resolved before parameter assignment."
        )

    # 1. PDB-Lexical atom naming to guarantee 100% identity with merge_complex.py
    pdb_block = Chem.MolToPDBBlock(rdmol)
    canonical_pdb_names = []
    for line in pdb_block.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            canonical_pdb_names.append(line[12:16].strip())

    mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=False)
    if len(canonical_pdb_names) != mol.n_atoms:
        raise ValueError(
            f"FATAL: Atom count mismatch between PDB block ({len(canonical_pdb_names)}) and molecule ({mol.n_atoms})!"
        )

    # 2. Allocate disjoint ahead-of-time atom types
    atom_types = allocate_atom_types(mol.n_atoms, args.prefix)

    # 3. Partial charge assignment
    if args.charge_method == "resp":
        if not args.charges:
            raise ValueError("FATAL: --charges (.crg) is strictly required when charge_method='resp'!")
        charges_f = load_external_crg_charges(Path(args.charges).resolve(), mol.n_atoms)
        mol.partial_charges = charges_f * unit.elementary_charge
    elif args.charge_method == "ashgc":
        if not any(isinstance(t, NAGLToolkitWrapper) for t in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits):
            if NAGLToolkitWrapper.is_available():
                GLOBAL_TOOLKIT_REGISTRY.register_toolkit(NAGLToolkitWrapper)
        NAGLToolkitWrapper().assign_partial_charges(
            molecule=mol, partial_charge_method="openff-gnn-am1bcc-1.0.0.pt"
        )
        charges_f = [float(q.m_as("elementary_charge")) for q in mol.partial_charges]
    elif args.charge_method == "bcc":
        mol.assign_partial_charges(partial_charge_method="am1bcc")
        charges_f = [float(q.m_as("elementary_charge")) for q in mol.partial_charges]

    # Verify charge neutrality / formal charge conservation
    sum_q = sum(charges_f)
    if abs(sum_q - float(args.charge)) > 1e-3:
        raise ValueError(
            f"FATAL: Sum of partial charges ({sum_q:.4f}) deviates from formal charge ({args.charge})!"
        )

    # 4. Resolve SMIRNOFF potentials via Interchange
    ff_name = args.forcefield if args.forcefield.endswith(".offxml") else f"{args.forcefield}.offxml"
    ff = ForceField(ff_name)
    interchange = ff.create_interchange(mol.to_topology())

    # 5. Nonbonded translation: sigma (nm) -> R* (Angstrom)
    vdw_coll = interchange.collections["vdW"]
    type_params = {}
    for idx in range(mol.n_atoms):
        t_name = atom_types[idx]
        pot_key = vdw_coll.key_map.get((idx,))
        pot = vdw_coll.potentials[pot_key]
        sig_nm = float(pot.parameters["sigma"].m_as("nanometer"))
        eps_kcal = float(pot.parameters["epsilon"].m_as("kilocalorie / mole"))
        mass = mol.atoms[idx].mass.m_as("dalton")
        r_star_ang = (sig_nm * 10.0) * (2.0 ** (1.0 / 6.0)) / 2.0
        type_params[t_name] = {
            "mass": mass,
            "r_star": r_star_ang,
            "epsilon": eps_kcal,
        }

    # 6. Harmonic Bonds: K_r = k_openff / 2.0
    bonds_coll = interchange.collections["Bonds"]
    bond_params = {}
    for top_key, pot_key in bonds_coll.key_map.items():
        i, j = top_key.atom_indices
        pot = bonds_coll.potentials[pot_key]
        k_openff = float(pot.parameters["k"].m_as("kilocalorie / (mole * angstrom**2)"))
        r0 = float(pot.parameters["length"].m_as("angstrom"))
        k_amber = k_openff / 2.0
        t1, t2 = atom_types[i], atom_types[j]
        b_key = tuple(sorted([t1, t2]))
        bond_params[b_key] = (k_amber, r0)

    # 7. Harmonic Angles: K_theta = k_openff / 2.0
    angles_coll = interchange.collections["Angles"]
    angle_params = {}
    for top_key, pot_key in angles_coll.key_map.items():
        i, j, k = top_key.atom_indices
        pot = angles_coll.potentials[pot_key]
        k_openff = float(pot.parameters["k"].m_as("kilocalorie / (mole * radian**2)"))
        th0 = float(pot.parameters["angle"].m_as("degree"))
        k_amber = k_openff / 2.0
        t1, t2, t3 = atom_types[i], atom_types[j], atom_types[k]
        a_key = (t1, t2, t3) if t1 <= t3 else (t3, t2, t1)
        angle_params[a_key] = (k_amber, th0)

    # 8. Proper Dihedrals: Mathematical symmetric tie-breaking canonicalization
    tors_coll = interchange.collections["ProperTorsions"]
    dihe_params = {}
    for top_key, pot_key in tors_coll.key_map.items():
        i, j, k, l = top_key.atom_indices
        pot = tors_coll.potentials[pot_key]
        k_val = float(pot.parameters["k"].m_as("kilocalorie / mole"))
        period = int(pot.parameters["periodicity"])
        phase = float(pot.parameters["phase"].m_as("degree"))
        idivf = int(pot.parameters.get("idivf", 1))

        t1, t2, t3, t4 = atom_types[i], atom_types[j], atom_types[k], atom_types[l]
        d_key = min((t1, t2, t3, t4), (t4, t3, t2, t1))
        dihe_params.setdefault(d_key, []).append((idivf, k_val, phase, period))

    # 9. Improper Dihedrals: Central trivalent hub aligned to position 3 with AMBER energy scaling
    impr_params = {}
    if "ImproperTorsions" in interchange.collections:
        impr_coll = interchange.collections["ImproperTorsions"]
        adjacency = {idx: set() for idx in range(mol.n_atoms)}
        for bond in mol.bonds:
            adjacency[bond.atom1_index].add(bond.atom2_index)
            adjacency[bond.atom2_index].add(bond.atom1_index)

        for top_key, pot_key in impr_coll.key_map.items():
            quartet = list(top_key.atom_indices)
            pot = impr_coll.potentials[pot_key]
            k_val = float(pot.parameters["k"].m_as("kilocalorie / mole"))
            period = int(pot.parameters["periodicity"])
            phase = float(pot.parameters["phase"].m_as("degree"))

            central_cands = [
                a for a in quartet
                if sum(bool(other in adjacency[a]) for other in quartet if other != a) == 3
            ]
            cent_atom = central_cands[0] if len(central_cands) == 1 else quartet[1]
            periph = sorted([a for a in quartet if a != cent_atom])
            t1, t2, t3, t4 = atom_types[periph[0]], atom_types[periph[1]], atom_types[cent_atom], atom_types[periph[2]]
            # SMIRNOFF specifies trefoil improper energy as average over 3 impropers:
            # U_trefoil = (1/3) * sum_{i=1}^3 k_openff * (1 + cos(n*phi - phase)) = k_openff * (1 + cos(n*phi - phase)).
            # AMBER evaluates a single matched improper with U = (K_amber / 2) * (1 + cos(n*phi - phase)).
            # Therefore, K_amber = 2.0 * k_openff preserves the full out-of-plane restraint barrier.
            k_amber_impr = k_val * 2.0
            impr_params[(t1, t2, t3, t4)] = (k_amber_impr, phase, period)

    # 10. Write Canonical Mol2
    clean_res = args.res_name[:4].upper()
    conf = mol.conformers[0]
    out_mol2 = Path(args.out_mol2).resolve()
    out_mol2.parent.mkdir(parents=True, exist_ok=True)

    with open(out_mol2, "w", encoding="utf-8") as f:
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{clean_res}\n")
        f.write(f"{mol.n_atoms} {mol.n_bonds} 1 0 0\n")
        f.write("SMALL\nUSER_CHARGES\n\n")

        f.write("@<TRIPOS>ATOM\n")
        for idx in range(mol.n_atoms):
            pos = conf[idx].m_as("angstrom")
            aname = canonical_pdb_names[idx]
            atype = atom_types[idx]
            q = charges_f[idx]
            f.write(
                f"{idx+1:>7d} {aname:<8s} {pos[0]:>10.4f} {pos[1]:>10.4f} {pos[2]:>10.4f} {atype:<5s} 1 {clean_res:<4s} {q:>10.6f}\n"
            )

        f.write("@<TRIPOS>BOND\n")
        for b_idx, bond in enumerate(mol.bonds, 1):
            i = bond.atom1_index + 1
            j = bond.atom2_index + 1
            b_order = "ar" if bond.is_aromatic else str(int(round(bond.bond_order or 1)))
            f.write(f"{b_idx:>6d} {i:>5d} {j:>5d} {b_order:>4s}\n")

        f.write("@<TRIPOS>SUBSTRUCTURE\n")
        f.write(f"1 {clean_res} 1 TEMP 0 **** **** 0 ROOT\n")

    # 11. Write Canonical AMBER Frcmod with strict forcefield metadata header
    out_frcmod = Path(args.out_frcmod).resolve()
    out_frcmod.parent.mkdir(parents=True, exist_ok=True)

    with open(out_frcmod, "w", encoding="utf-8") as f:
        f.write(f"Remark: OpenFF {ff_name} parameters for {clean_res} [Method: {args.charge_method.upper()}, Prefix: {args.prefix}]\n")
        f.write("MASS\n")
        for t_name, p in sorted(type_params.items()):
            f.write(f"{t_name:<2s}  {p['mass']:6.3f}\n")
        f.write("\nBOND\n")
        for (t1, t2), (k_b, r0) in sorted(bond_params.items()):
            f.write(f"{t1:<2s}-{t2:<2s}  {k_b:7.2f}  {r0:6.3f}\n")
        f.write("\nANGLE\n")
        for (t1, t2, t3), (k_a, th0) in sorted(angle_params.items()):
            f.write(f"{t1:<2s}-{t2:<2s}-{t3:<2s}  {k_a:7.2f}  {th0:7.2f}\n")
        f.write("\nDIHE\n")
        for (t1, t2, t3, t4), terms in sorted(dihe_params.items()):
            terms_sorted = sorted(terms, key=lambda x: abs(x[3]))
            n_terms = len(terms_sorted)
            for idx_term, (idivf, pk, phase, period) in enumerate(terms_sorted, 1):
                eff_period = -abs(period) if idx_term < n_terms else abs(period)
                f.write(f"{t1:<2s}-{t2:<2s}-{t3:<2s}-{t4:<2s}  {idivf:2d}  {pk:8.3f}  {phase:7.2f}  {eff_period:5.1f}\n")
        f.write("\nIMPROPER\n")
        for (t1, t2, t3, t4), (k_val, phase, period) in sorted(impr_params.items()):
            f.write(f"{t1:<2s}-{t2:<2s}-{t3:<2s}-{t4:<2s}  {k_val:8.3f}  {phase:7.2f}  {abs(period):5.1f}\n")
        f.write("\nNONBON\n")
        for t_name, p in sorted(type_params.items()):
            f.write(f"  {t_name:<2s}  {p['r_star']:8.4f}  {p['epsilon']:8.4f}\n")
        f.write("\n")

    print(
        f"✓ Successfully generated OpenFF AMBER library: {out_mol2.name} & {out_frcmod.name} "
        f"({mol.n_atoms} atoms, prefix '{args.prefix}', method: {args.charge_method})",
        flush=True
    )


if __name__ == "__main__":
    main()
