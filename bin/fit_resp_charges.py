#!/usr/bin/env python3
"""
bin/fit_resp_charges.py
Extracts Merz-Kollman ESP charges from ORCA output, enforces topological symmetry
across graph automorphisms via RDKit, and guarantees exact net charge conservation.
"""

import argparse
import sys
from pathlib import Path
from rdkit import Chem


def parse_orca_resp_charges(orca_out_path: Path, expected_atoms: int) -> list[float]:
    """Parses raw Merz-Kollman RESP charges from ORCA standard output."""
    in_resp_section = False
    charges_dict = {}

    with open(orca_out_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if "RESP CHARGES" in line_str:
                in_resp_section = True
                continue
            if in_resp_section:
                if (
                    line_str.startswith("ATOM")
                    or line_str.startswith("-")
                    or line_str.startswith("=")
                ):
                    continue
                if (
                    not line_str
                    or line_str.startswith("Sum of")
                    or line_str.startswith("---")
                ):
                    if charges_dict:
                        break
                    continue
                clean_line = line_str.replace(":", " ")
                parts = clean_line.split()
                if len(parts) >= 3 and parts[0].isdigit():
                    idx = int(parts[0])
                    val = float(parts[2])
                    charges_dict[idx] = val

    if len(charges_dict) != expected_atoms:
        raise ValueError(
            f"FATAL: Expected {expected_atoms} atomic charges from ORCA RESP output, "
            f"but parsed {len(charges_dict)} from {orca_out_path}!"
        )

    return [charges_dict[i] for i in range(expected_atoms)]


def symmetrize_charges(
    mol: Chem.Mol, raw_charges: list[float], formal_charge: int
) -> list[float]:
    """
    Symmetrizes partial atomic charges across topological automorphisms
    and guarantees exact integer net charge conservation.
    """
    n_atoms = mol.GetNumAtoms()
    matches = mol.GetSubstructMatches(mol, uniquify=False)

    equiv_classes = []
    visited = set()
    for i in range(n_atoms):
        if i in visited:
            continue
        eq_group = {match[i] for match in matches}
        visited.update(eq_group)
        equiv_classes.append(sorted(list(eq_group)))

    symmetrized = list(raw_charges)
    for group in equiv_classes:
        mean_q = sum(raw_charges[idx] for idx in group) / len(group)
        for idx in group:
            symmetrized[idx] = mean_q

    # Enforce exact net charge conservation without rounding drift
    total_q = sum(symmetrized)
    drift = float(formal_charge) - total_q
    correction_per_atom = drift / float(n_atoms)

    final_charges = [q + correction_per_atom for q in symmetrized]
    final_sum = sum(final_charges)

    if abs(final_sum - float(formal_charge)) > 1e-5:
        raise ValueError(
            f"FATAL: Symmetrized charge sum ({final_sum:.6f}) deviates from "
            f"formal charge ({formal_charge}) beyond physical tolerance!"
        )

    return final_charges


def main():
    parser = argparse.ArgumentParser(
        description="Fit and symmetrize atomic ESP charges"
    )
    parser.add_argument("--sdf", required=True, help="Input relaxed 3D structure")
    parser.add_argument(
        "--orca_log", required=True, help="ORCA ESP calculation output log"
    )
    parser.add_argument("--charge", required=True, type=int, help="Net formal charge")
    parser.add_argument(
        "--out_crg", required=True, help="Output vector file of atomic charges (.crg)"
    )
    parser.add_argument("--out_sdf", required=True, help="Output charged SDF structure")
    args = parser.parse_args()

    sdf_path = Path(args.sdf).resolve()
    orca_log = Path(args.orca_log).resolve()

    if not sdf_path.exists():
        raise FileNotFoundError(f"FATAL: SDF structure not found: {sdf_path}")
    if not orca_log.exists():
        raise FileNotFoundError(f"FATAL: ORCA log not found: {orca_log}")

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(
            f"FATAL: Could not read valid molecular structure from {sdf_path}"
        )

    n_atoms = mol.GetNumAtoms()
    raw_charges = parse_orca_resp_charges(orca_log, n_atoms)
    final_charges = symmetrize_charges(mol, raw_charges, args.charge)

    # Write clean flat charges file formatted for antechamber -cf
    out_crg_path = Path(args.out_crg).resolve()
    out_crg_path.parent.mkdir(parents=True, exist_ok=True)
    out_crg_path.write_text(
        "\n".join(f"{q:14.6f}" for q in final_charges) + "\n", encoding="utf-8"
    )

    out_sdf_path = Path(args.out_sdf).resolve()
    out_sdf_path.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(out_sdf_path)) as writer:
        writer.SetForceV3000(False)
        writer.write(mol)

    print(
        f"✓ Successfully extracted and symmetrized {n_atoms} RESP charges -> {out_crg_path.name}",
        flush=True,
    )


if __name__ == "__main__":
    main()
