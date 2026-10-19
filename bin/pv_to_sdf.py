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
            content = f.read()
        
        blocks = content.split('$$$$\n')
        with open(args.sdf, "w", encoding="utf-8") as f:
            for block in blocks:
                if not block.strip():
                    continue
                lines = block.split('\n')
                # Line 2 is the program/timestamp/dimension line in MDL format
                if len(lines) > 1:
                    if len(lines[1]) >= 22 and lines[1][20:22] == '2D':
                        lines[1] = lines[1][:20] + '3D' + lines[1][22:]
                    elif '2D' in lines[1]:
                        lines[1] = lines[1].replace('2D', '3D')
                f.write('\n'.join(lines) + '$$$$\n')
    except Exception as e:
        import sys
        print(f"Warning: Could not format dimensionality tags in {args.sdf}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
