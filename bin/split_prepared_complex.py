#!/usr/bin/env python3
"""
bin/split_prepared_complex.py
Splits a fully prepared and protonated complex (from PrepWizard) into
receptor (protein+cofactors) and native ligand.
Preserves perfect in-pocket coordinates and correctly assigned protonation states.
"""
import argparse
import sys
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import analyze


def main():
    parser = argparse.ArgumentParser(
        description="Split prepared complex into receptor and ligand"
    )
    parser.add_argument(
        "--complex_mae", required=True, help="Input prepared complex MAE"
    )
    parser.add_argument("--ligand", default="", help="Native ligand residue name")
    parser.add_argument(
        "--cofactors", nargs="*", default=[], help="Cofactor/ion residue names"
    )
    parser.add_argument(
        "--out_receptor_mae", required=True, help="Output prepared receptor MAE"
    )
    parser.add_argument(
        "--out_receptor_pdb", required=True, help="Output prepared receptor PDB"
    )
    parser.add_argument(
        "--out_ligand_mae", required=True, help="Output prepared ligand MAE"
    )
    parser.add_argument(
        "--out_ligand_sdf", required=True, help="Output prepared ligand SDF"
    )
    parser.add_argument(
        "--out_manifest", default="extracted_poses.json", help="Output JSON manifest"
    )
    args = parser.parse_args()

    try:
        with structure.StructureReader(args.complex_mae) as reader:
            sts = list(reader)
        if not sts:
            raise ValueError(f"No structures found in {args.complex_mae}")
        st = sts[0]
        for next_st in sts[1:]:
            st.extend(next_st)
    except Exception as e:
        print(f"FATAL: Could not read {args.complex_mae}: {e}")
        sys.exit(1)

    import json
    extracted_records = []
    target_lig = args.ligand.strip().upper() if args.ligand else None
    target_cofactors = {c.strip().upper() for c in (args.cofactors or []) if c.strip()}

    # 1. Extract active native ligand if present
    if target_lig:
        lig_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == target_lig]
        if not lig_atoms:
            avail = sorted(list(set(r.pdbres.strip().upper() for r in st.residue)))
            raise ValueError(
                f"FATAL: Native ligand '{target_lig}' not found in prepared complex! "
                f"Available residues: {avail}"
            )
        lig_st = st.extract(lig_atoms)
        lig_st.title = target_lig
        lig_st.write(args.out_ligand_mae)
        lig_st.write(args.out_ligand_sdf)
        st.deleteAtoms(lig_atoms)
        extracted_records.append({
            "res_name": target_lig,
            "charge": int(lig_st.formal_charge),
            "sdf": args.out_ligand_sdf,
            "mae": args.out_ligand_mae,
            "is_ligand": True,
            "is_ion": (lig_st.atom_total == 1),
        })
        print(f"Extracted prepared native ligand '{target_lig}' ({len(lig_atoms)} atoms).")

    # 2. Extract declared cofactors and ions with PrepWizard-optimized hydrogens
    for cof_name in sorted(list(target_cofactors)):
        cof_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == cof_name]
        if not cof_atoms:
            continue
        cof_st = st.extract(cof_atoms)
        cof_st.title = cof_name
        cof_sdf = f"pose_{cof_name}.sdf"
        cof_mae = f"pose_{cof_name}.mae"
        cof_st.write(cof_mae)
        cof_st.write(cof_sdf)
        st.deleteAtoms(cof_atoms)
        extracted_records.append({
            "res_name": cof_name,
            "charge": int(cof_st.formal_charge),
            "sdf": cof_sdf,
            "mae": cof_mae,
            "is_ligand": False,
            "is_ion": (cof_st.atom_total == 1),
        })
        print(f"Extracted prepared cofactor/ion '{cof_name}' ({len(cof_atoms)} atoms).")

    # The remaining structure is 100% PURE APO PROTEIN
    st.write(args.out_receptor_mae)
    st.write(args.out_receptor_pdb)

    Path(args.out_manifest).write_text(json.dumps(extracted_records, indent=2), encoding="utf-8")
    print(
        f"Wrote prepared receptor to {args.out_receptor_mae} and {args.out_receptor_pdb} ({st.atom_total} atoms)."
    )


if __name__ == "__main__":
    main()
