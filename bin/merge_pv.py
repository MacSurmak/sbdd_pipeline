#!/usr/bin/env python3
"""
bin/merge_pv.py
Safely merges chunked Glide Pose Viewer (_pv.maegz) files into a single valid PV file,
preserving the receptor strictly once at index 0 and appending all ligand poses.
Also merges corresponding CSV and SDF pose files.
"""
import argparse
import csv
import sys
from pathlib import Path
from schrodinger import structure


def main():
    parser = argparse.ArgumentParser(
        description="Merge Glide Pose Viewer chunks without receptor duplication"
    )
    parser.add_argument(
        "--pv_files", nargs="+", required=True, help="Input chunk _pv.maegz files"
    )
    parser.add_argument(
        "--csv_files", nargs="*", default=[], help="Input chunk score CSV files"
    )
    parser.add_argument(
        "--sdf_files", nargs="*", default=[], help="Input chunk pose SDF files"
    )
    parser.add_argument("--out_pv", required=True, help="Output merged _pv.maegz path")
    parser.add_argument("--out_csv", help="Output merged CSV path")
    parser.add_argument("--out_sdf", help="Output merged SDF path")
    args = parser.parse_args()

    receptor_written = False
    total_ligand_poses = 0
    pv_files = sorted(args.pv_files)

    # 1. Merge PV files
    with structure.StructureWriter(args.out_pv) as writer:
        for pv_file in pv_files:
            pv_path = Path(pv_file).resolve()
            if not pv_path.exists() or pv_path.stat().st_size == 0:
                print(f"FATAL: PV chunk file is missing or empty: {pv_path}", file=sys.stderr)
                sys.exit(1)

            reader = None
            try:
                reader = structure.StructureReader(str(pv_path))
                for idx, st in enumerate(reader):
                    is_receptor = (idx == 0) and (
                        st.property.get("b_glide_receptor", 0) == 1
                        or st.atom_total > 500
                    )
                    if is_receptor:
                        if not receptor_written:
                            writer.append(st)
                            receptor_written = True
                        continue

                    # Append ligand pose
                    writer.append(st)
                    total_ligand_poses += 1
            finally:
                if reader is not None:
                    reader.close()

    print(
        f"Merged {len(pv_files)} PV files: 1 receptor + {total_ligand_poses} ligand poses -> {args.out_pv}"
    )

    # 2. Merge SDF files if requested
    if args.out_sdf and args.sdf_files:
        with open(args.out_sdf, "w", encoding="utf-8") as out_sdf_f:
            for sdf_file in sorted(args.sdf_files):
                p = Path(sdf_file)
                if not p.exists() or p.stat().st_size == 0:
                    print(f"FATAL: SDF chunk file is missing or empty: {sdf_file}", file=sys.stderr)
                    sys.exit(1)
                out_sdf_f.write(p.read_text(encoding="utf-8", errors="ignore"))
        print(f"Merged SDF pose files -> {args.out_sdf}")

    # 3. Merge CSV score files if requested
    if args.out_csv and args.csv_files:
        header = None
        rows = []
        for csv_file in sorted(args.csv_files):
            p = Path(csv_file)
            if not p.exists() or p.stat().st_size == 0:
                print(f"FATAL: CSV chunk file is missing or empty: {csv_file}", file=sys.stderr)
                sys.exit(1)
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                r = csv.reader(f)
                for r_idx, row in enumerate(r):
                    if r_idx == 0:
                        if header is None:
                            header = row
                    else:
                        rows.append(row)

        with open(args.out_csv, "w", encoding="utf-8", newline="") as out_csv_f:
            w = csv.writer(out_csv_f)
            if header:
                w.writerow(header)
            w.writerows(rows)
        print(f"Merged {len(rows)} score records -> {args.out_csv}")


if __name__ == "__main__":
    main()
