#!/usr/bin/env python3
"""
bin/select_hits.py
Extracts top-ranked hits from ranking.csv using standard library csv.
"""
import argparse
import csv
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Select top hits from ranking.csv")
    parser.add_argument("--ranking_csv", required=True, help="Input ranking.csv file")
    parser.add_argument(
        "--top_n", type=int, default=1, help="Number of top hits to select"
    )
    parser.add_argument(
        "--include_native",
        default="true",
        help="Always include native reference ligand",
    )
    parser.add_argument(
        "--native_ligand", default="", help="Exact title of the native ligand"
    )
    parser.add_argument(
        "--out",
        default="selected_hits.txt",
        help="Output text file with selected titles",
    )
    args = parser.parse_args()

    selected = []
    csv_path = Path(args.ranking_csv)
    if csv_path.exists() and csv_path.stat().st_size > 0:
        with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f)
            for row in reader:
                title = (
                    row.get("Title") or row.get("title") or row.get("s_m_title") or ""
                )
                title = title.strip()
                frame = row.get("Frame") or row.get("frame")

                if not frame:
                    raise ValueError(
                        f"FATAL: Missing 'Frame' column in ranking row for {title}"
                    )

                frame_id = str(frame).strip()
                if not frame_id:
                    raise ValueError(
                        f"FATAL: Empty 'Frame' value in ranking row for {title}"
                    )

                if title and not any(s[0] == title for s in selected):
                    selected.append((title, frame_id))
                    if len(selected) >= args.top_n:
                        break
    else:
        raise FileNotFoundError(
            f"FATAL: Ranking CSV {args.ranking_csv} not found or empty."
        )

    include_nat = str(args.include_native).lower() in ["true", "1", "yes"]
    nat_lig = args.native_ligand.strip()

    if include_nat and nat_lig and not any(s[0] == nat_lig for s in selected):
        selected.append((nat_lig, "1"))

    if not selected:
        raise ValueError("FATAL: No hits selected and no native ligand specified.")

    print(f"Selected {len(selected)} hit(s) for MD simulation: {selected}")
    with open(args.out, "w", encoding="utf-8") as f:
        for title, fid in selected:
            f.write(f"{title},{fid}\n")


if __name__ == "__main__":
    main()
