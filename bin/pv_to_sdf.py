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

    # Physically correct Schrodinger's invalid 2D tags in the MDL V2000 header to strictly 
    # prevent RDKit coordinate dimensionality warnings without silencing the logger.
    try:
        with open(args.sdf, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
        for i in range(len(lines)):
            # Line 1 (0-indexed) is the timestamp/dimension line. For subsequent molecules, it's 2 lines after $$$$
            if i == 1 or (i > 1 and lines[i-2] == '$$$$'):
                if len(lines[i]) >= 22 and lines[i][20:22] == '2D':
                    lines[i] = lines[i][:20] + '3D' + lines[i][22:]
                elif '2D' in lines[i]:
                    lines[i] = lines[i].replace('2D', '3D')
        with open(args.sdf, "w", encoding="utf-8") as f:
            f.write('\n'.join(lines) + '\n')
    except Exception as e:
        import sys
        print(f"Warning: Could not format dimensionality tags in {args.sdf}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
