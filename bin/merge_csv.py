#!/usr/bin/env python3
"""
bin/merge_csv.py
Safely merges multiple CSV files with identical headers.
"""
import argparse
import csv
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True, help="Input CSV files")
    parser.add_argument("--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    header = None
    rows = []

    for f in sorted(args.inputs):
        p = Path(f)
        if not p.exists() or p.stat().st_size == 0:
            print(f"FATAL: CSV file is missing or empty: {p}", file=sys.stderr)
            sys.exit(1)

        with open(p, "r", encoding="utf-8", errors="ignore") as fp:
            r = csv.reader(fp)
            for i, row in enumerate(r):
                if i == 0:
                    if header is None:
                        header = row
                else:
                    rows.append(row)

    with open(args.output, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out)
        if header:
            w.writerow(header)
        w.writerows(rows)


if __name__ == "__main__":
    main()
