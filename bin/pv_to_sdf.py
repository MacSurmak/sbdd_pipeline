"""
Extracts docked ligand poses from Pose Viewer MAE to a multi-conformation SDF.
"""

import argparse
from schrodinger import structure


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pv", required=True, help="Input _pv.maegz file")
    parser.add_argument("--sdf", required=True, help="Output .sdf file")
    parser.add_argument(
        "--max_poses", type=int, default=10, help="Maximum poses per ligand to export"
    )
    args = parser.parse_args()

    pose_counts = {}
    with (
        structure.StructureReader(args.pv) as reader,
        structure.StructureWriter(args.sdf) as writer,
    ):
        for i, st in enumerate(reader):
            # Skip the receptor
            if i == 0 and (
                st.property.get("b_glide_receptor", 0) == 1 or st.atom_total > 500
            ):
                continue
            title = (st.property.get("s_m_title") or st.title or "").strip()
            # Retain only up to max_poses for each unique ligand title
            if not title or pose_counts.get(title, 0) < args.max_poses:
                writer.append(st)
                if title:
                    pose_counts[title] = pose_counts.get(title, 0) + 1

    # Structurally normalize the SDF via RDKit to ensure proper 3D MDL V2000 headers
    # and eliminate downstream dimensionality warnings natively. Fail-Fast on error.
    try:
        import os
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")

        suppl = Chem.SDMolSupplier(args.sdf, removeHs=False, sanitize=False)
        tmp_sdf = args.sdf + ".tmp"
        with Chem.SDWriter(tmp_sdf) as w:
            for m in suppl:
                if m is not None:
                    w.write(m)
        os.replace(tmp_sdf, args.sdf)
    except Exception as e:
        import sys

        print(
            f"FATAL: Failed to structurally normalize SDF {args.sdf} via RDKit: {e}",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
