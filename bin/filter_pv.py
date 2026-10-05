"""
Worker script executed via $SCHRODINGER/run python3.
Filters a Pose Viewer file, keeping only ligands explicitly present in the whitelist for a given frame.
Supports splitting the filtered poses into chunks of N poses, each with the receptor prepended,
enabling massive distributed Prime MM-GBSA across cluster nodes.
"""

import argparse
import sys
import csv
from schrodinger import structure


def main():
    parser = argparse.ArgumentParser(
        description="Filter Pose Viewer file by whitelist and optionally chunk"
    )
    parser.add_argument("--input", required=True, help="Input _pv.maegz file")
    parser.add_argument("--whitelist", required=True, help="Path to whitelist.csv")
    parser.add_argument("--frame", required=True, help="Current frame ID")
    parser.add_argument(
        "--max_poses", type=int, default=1, help="Maximum poses per ligand for MM-GBSA"
    )
    parser.add_argument(
        "--out_prefix",
        default="",
        help="Prefix for chunked outputs, e.g. frame_filtered",
    )
    parser.add_argument(
        "--batch_size",
        "--chunk_size",
        dest="batch_size",
        type=int,
        default=25,
        help="Poses per chunk",
    )
    parser.add_argument(
        "--output", default="", help="Optional single output _pv.maegz file"
    )
    args = parser.parse_args()

    import re

    allowed_titles = set()
    with open(args.whitelist, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            f_clean = str(row["Frame"]).replace("F", "").replace("f", "").strip()
            req_clean = str(args.frame).replace("F", "").replace("f", "").strip()
            if f_clean == req_clean:
                raw_t = row["Title"].strip()
                allowed_titles.add(raw_t)
                allowed_titles.add(re.sub(r"(_model_\d+|_dock_\d+|_\d+)$", "", raw_t))

    receptor = None
    whitelisted_poses = []
    pose_counts = {}

    with structure.StructureReader(args.input) as reader:
        for i, st in enumerate(reader):
            if i == 0:
                receptor = st.copy()
                continue

            title = (st.property.get("s_m_title") or st.title or "").strip()
            title_norm = re.sub(r"(_model_\d+|_dock_\d+|_\d+)$", "", title)
            is_match = (title in allowed_titles) or (title_norm in allowed_titles)
            if is_match and pose_counts.get(title, 0) < args.max_poses:
                whitelisted_poses.append(st.copy())
                pose_counts[title] = pose_counts.get(title, 0) + 1

    if not whitelisted_poses:
        print(
            f"FATAL: Frame {args.frame}: No whitelisted poses found. Filter returned empty result.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(
        f"Frame {args.frame}: Identified receptor and {len(whitelisted_poses)} whitelisted poses."
    )

    # 1. Single file output if requested
    if args.output:
        with structure.StructureWriter(args.output) as writer:
            if receptor:
                writer.append(receptor)
            for p in whitelisted_poses:
                writer.append(p)
        print(f"Wrote single filtered PV to {args.output}")

    # Ensure receptor is explicitly flagged as glide receptor for Prime
    if receptor:
        receptor.property["b_glide_receptor"] = 1

    # 2. Chunked output if out_prefix specified
    if args.out_prefix:
        if not whitelisted_poses:
            print(
                f"Frame {args.frame}: No whitelisted poses found. Skipping chunk generation."
            )
        else:
            chunk_idx = 0
            for i in range(0, len(whitelisted_poses), args.batch_size):
                chunk = whitelisted_poses[i : i + args.batch_size]
                chunk_file = f"{args.out_prefix}_c{chunk_idx}.maegz"
                with structure.StructureWriter(chunk_file) as writer:
                    if receptor:
                        writer.append(receptor)
                    for p in chunk:
                        writer.append(p)
                print(
                    f"Wrote chunk {chunk_idx}: receptor + {len(chunk)} poses -> {chunk_file}"
                )
                chunk_idx += 1


if __name__ == "__main__":
    main()
