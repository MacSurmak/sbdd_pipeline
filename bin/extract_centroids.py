#!/usr/bin/env python3
"""
bin/extract_centroids.py
Atomic Stage 3: Direct seek frame extraction of physical medoid structures.
Preserves true chemical graph topology (formal charges, bond orders, stereochemistry)
via coordinate transfer onto pristine SubstanceState templates.
Extracts pure receptor, complex view, and all individual cofactor/ion poses.
"""

import argparse
import json
import sys
from pathlib import Path
import MDAnalysis as mda
from rdkit import Chem


STANDARD_PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLN", "GLU", "GLY", "HIS",
    "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL", "ASH", "GLH", "ACE", "NME", "NMA",
}

WATER_RESIDUES = {"HOH", "WAT", "SOL", "TIP", "TIP3", "TIP4", "W"}


def transfer_coordinates_to_template(template_sdf: Path, raw_mol: Chem.Mol, out_sdf: Path) -> None:
    """
    Transfers exact Cartesian coordinates from MD frame onto canonical SubstanceState template
    using strict topological graph isomorphism matching.
    Guarantees that atom serials in PDB do not scramble coordinates across different atom types.
    """
    suppl = Chem.SDMolSupplier(str(template_sdf), removeHs=False, sanitize=False)
    tmpl_mol = next((m for m in suppl if m is not None), None)
    if not tmpl_mol or not raw_mol:
        raise ValueError(f"FATAL: Failed to parse template {template_sdf} or raw MD conformer!")

    n_tmpl = tmpl_mol.GetNumAtoms()
    n_raw = raw_mol.GetNumAtoms()
    if n_tmpl != n_raw:
        raise ValueError(
            f"FATAL: Atom count mismatch between template ({n_tmpl}) and raw pose ({n_raw}) for {out_sdf.name}!"
        )

    # Build pure topology queries to find unambiguous isomorphism mapping
    q_tmpl = Chem.Mol(tmpl_mol)
    q_raw = Chem.Mol(raw_mol)
    for q in (q_tmpl, q_raw):
        for b in q.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        for a in q.GetAtoms():
            a.SetFormalCharge(0)
            a.SetIsotope(0)

    matches = q_raw.GetSubstructMatches(q_tmpl, uniquify=False)
    if not matches:
        raise ValueError(
            f"FATAL: Topological graph isomorphism failed! MD frame molecule in {out_sdf.name} "
            f"does not match template graph {template_sdf.name}."
        )

    match = matches[0]
    conf_tmpl = tmpl_mol.GetConformer()
    conf_raw = raw_mol.GetConformer()

    for tmpl_idx, raw_idx in enumerate(match):
        pos = conf_raw.GetAtomPosition(raw_idx)
        if abs(pos.x) < 1e-4 and abs(pos.y) < 1e-4 and abs(pos.z) < 1e-4:
            raise ValueError(
                f"FATAL: Degenerate coordinate (0, 0, 0) detected at atom index {raw_idx} during transfer!"
            )
        conf_tmpl.SetAtomPosition(tmpl_idx, pos)

    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(out_sdf)) as writer:
        writer.SetForceV3000(False)
        writer.write(tmpl_mol)


def main():
    parser = argparse.ArgumentParser(description="Extract validated 3D medoid structures from trajectory")
    parser.add_argument("--struct", required=True, help="Clean reference structure PDB")
    parser.add_argument("--xtc", nargs="+", required=True, help="Clean trajectory XTC file(s)")
    parser.add_argument("--medoids_json", required=True, help="Input medoids_summary.json deliverable")
    parser.add_argument("--substance_template_sdf", default="", help="Optional active ligand SubstanceState template")
    parser.add_argument("--ligand_resname", default="", help="Active ligand residue name")
    parser.add_argument("--outdir", default="centroids_out", help="Output directory for structures")
    args = parser.parse_args()

    struct_p = Path(args.struct).resolve()
    medoids_p = Path(args.medoids_json).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not struct_p.exists() or not medoids_p.exists():
        raise FileNotFoundError("FATAL: Required reference structure or medoids summary missing!")

    clean_xtcs = []
    for raw_arg in args.xtc:
        for chunk in str(raw_arg).split():
            c_clean = chunk.strip()
            if c_clean:
                p_xtc = Path(c_clean).resolve()
                if p_xtc.exists() and p_xtc.stat().st_size > 0:
                    clean_xtcs.append(str(p_xtc))

    if not clean_xtcs:
        raise ValueError("FATAL: No valid trajectory files provided to extract_centroids!")

    with open(medoids_p, "r", encoding="utf-8") as f:
        medoids = json.load(f)

    u = mda.Universe(str(struct_p), clean_xtcs)
    protein_atoms = u.select_atoms("protein")
    if len(protein_atoms) == 0:
        raise ValueError("FATAL: No protein atoms found in reference structure!")

    # Discover all small molecules (active ligand, cofactors, ions) to preserve system integrity
    small_mol_resnames = set()
    for res in u.residues:
        r_name = res.resname.strip().upper()
        if r_name not in STANDARD_PROTEIN_RESIDUES and r_name not in WATER_RESIDUES:
            small_mol_resnames.add(r_name)

    tmpl_sdf_p = Path(args.substance_template_sdf).resolve() if args.substance_template_sdf else None
    active_res = args.ligand_resname.strip().upper() if args.ligand_resname else ""

    manifest = []
    for stat in medoids:
        rank = stat["rank"]
        m_frame = stat["medoid_frame"]
        pop_pct = stat["population_pct"]

        # Direct O(1) seek to frame
        u.trajectory[m_frame]

        # 1. Export pure apo protein receptor
        rec_pdb = outdir / f"centroid_{rank}_receptor.pdb"
        protein_atoms.write(str(rec_pdb))

        # 2. Export full assembled complex view for visualization
        view_pdb = outdir / f"centroid_{rank}_complex.pdb"
        u.atoms.write(str(view_pdb))

        extracted_poses = {}
        for r_name in sorted(list(small_mol_resnames)):
            mol_atoms = u.select_atoms(f"resname {r_name}")
            if len(mol_atoms) == 0:
                continue

            raw_pose_pdb = outdir / f".tmp_{rank}_{r_name}.pdb"
            mol_atoms.write(str(raw_pose_pdb))

            dest_sdf = outdir / f"centroid_{rank}_{r_name}_pose.sdf"
            raw_mol = Chem.MolFromPDBFile(str(raw_pose_pdb), removeHs=False, sanitize=False)
            if not raw_mol:
                raise ValueError(
                    f"FATAL: RDKit failed to parse extracted small molecule '{r_name}' from {raw_pose_pdb.name}!"
                )

            is_active_ligand = (r_name == active_res)
            if is_active_ligand and tmpl_sdf_p and tmpl_sdf_p.exists():
                transfer_coordinates_to_template(tmpl_sdf_p, raw_mol, dest_sdf)
            else:
                # Retain coordinates for monoatomic ions or cofactors
                with Chem.SDWriter(str(dest_sdf)) as w:
                    w.SetForceV3000(False)
                    w.write(raw_mol)

            raw_pose_pdb.unlink(missing_ok=True)
            extracted_poses[r_name] = {
                "res_name": r_name,
                "is_ligand": is_active_ligand,
                "sdf_path": str(dest_sdf.resolve()),
            }

        manifest.append({
            "rank": rank,
            "population_pct": pop_pct,
            "medoid_frame": m_frame,
            "medoid_time_ns": stat["medoid_time_ns"],
            "medoid_replicate": stat["medoid_replicate"],
            "receptor_pdb": str(rec_pdb.resolve()),
            "complex_pdb": str(view_pdb.resolve()),
            "poses": extracted_poses,
        })
        print(f"[EXTRACT_CENTROIDS] Centroid #{rank} extracted (frame {m_frame}, {pop_pct}%).", flush=True)

    with open(outdir / "extracted_manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"✓ Successfully extracted {len(manifest)} validated centroid structure deliverables.", flush=True)


if __name__ == "__main__":
    main()
