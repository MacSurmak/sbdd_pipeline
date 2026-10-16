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
        description="Deterministically split prepared complex into receptor and poses using manifest specification"
    )
    parser.add_argument(
        "--complex_mae", required=True, help="Input prepared complex MAE"
    )
    parser.add_argument(
        "--poses_spec_json",
        required=True,
        help="JSON file containing explicit pose extraction targets from DB plan",
    )
    parser.add_argument(
        "--out_receptor_mae",
        required=True,
        help="Explicit destination for prepared receptor MAE",
    )
    parser.add_argument(
        "--out_receptor_pdb",
        required=True,
        help="Explicit destination for prepared receptor PDB",
    )
    parser.add_argument(
        "--out_manifest",
        default="extracted_poses.json",
        help="Output JSON receipt of extracted poses",
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
        print(f"FATAL: Could not read {args.complex_mae}: {e}", file=sys.stderr)
        sys.exit(1)

    import json

    spec_path = Path(args.poses_spec_json)
    if not spec_path.exists():
        raise FileNotFoundError(
            f"FATAL: Pose specification manifest not found: {spec_path}"
        )

    poses_spec = json.loads(spec_path.read_text(encoding="utf-8"))
    extracted_records = []

    # Deterministic extraction strictly guided by DB specification
    for item in poses_spec:
        res_clean = item["res_name"].strip().upper()
        pose_name = item["pose_name"]
        inchikey = item["inchikey"]
        is_lig = item.get("is_ligand", False)
        is_ion = item.get("is_ion", False)

        atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == res_clean]
        if not atoms:
            avail = sorted(list(set(r.pdbres.strip().upper() for r in st.residue)))
            raise ValueError(
                f"FATAL: Residue '{res_clean}' for pose '{pose_name}' not found in prepared complex! "
                f"Available residues: {avail}"
            )

        pose_st = st.extract(atoms)

        # Strict Fail-Fast: Ensure organic poses are single connected graphs
        if not is_ion and len(pose_st.molecule) > 1:
            raise ValueError(
                f"FATAL: Residue '{res_clean}' resolves to {len(pose_st.molecule)} disconnected molecules! "
                f"Quantum parametrization (AM1-BCC) requires a single connected graph. Assign unique residue names."
            )

        pose_st.title = res_clean

        # Deterministic local filenames for Nextflow staging
        out_sdf = f"{pose_name}_pose.sdf"
        out_mae = f"{pose_name}_pose.mae"
        pose_st.write(out_sdf)
        pose_st.write(out_mae)

        st.deleteAtoms(atoms)

        extracted_records.append(
            {
                "pose_name": pose_name,
                "res_name": res_clean,
                "inchikey": inchikey,
                "charge": int(pose_st.formal_charge),
                "sdf": out_sdf,
                "mae": out_mae,
                "is_ligand": is_lig,
                "is_ion": is_ion,
            }
        )
        print(
            f"Extracted pose '{pose_name}' ({res_clean}, {len(atoms)} atoms, formal charge {int(pose_st.formal_charge)})."
        )

    # The remaining structure is 100% PURE APO PROTEIN
    st.write(args.out_receptor_mae)
    st.write(args.out_receptor_pdb)

    Path(args.out_manifest).write_text(
        json.dumps(extracted_records, indent=2), encoding="utf-8"
    )
    print(
        f"Wrote prepared receptor to {args.out_receptor_mae} and {args.out_receptor_pdb} ({st.atom_total} atoms)."
    )


if __name__ == "__main__":
    main()
