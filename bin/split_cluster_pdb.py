#!/usr/bin/env python3
"""
bin/split_cluster_pdb.py
Splits GROMACS multi-model cluster PDB into centroid structures and in-situ ligand poses.
Extracts both PDB and MAE for centroids, and SDF and MAE for ligand poses.
Outputs receipt TSV for CLI DB synchronization.
"""

import argparse
import sys
from pathlib import Path
from schrodinger import structure


def main():
    parser = argparse.ArgumentParser(
        description="Split cluster PDB and extract ligand poses"
    )
    parser.add_argument("--input", required=True, help="Input clusters_multi.pdb")
    parser.add_argument("--prefix", required=True, help="Prefix for output structures")
    parser.add_argument("--ligand", default="", help="Ligand residue name (e.g. KVR)")
    parser.add_argument(
        "--max_clusters", type=int, default=0, help="Max clusters to extract"
    )
    parser.add_argument(
        "--out_receipt",
        default="cluster_receipt_chunk.tsv",
        help="Output receipt TSV",
    )
    args = parser.parse_args()

    in_path = Path(args.input)
    if not in_path.exists() or in_path.stat().st_size == 0:
        raise FileNotFoundError(f"Cluster PDB {in_path} not found or empty.")

    lig_name = args.ligand.strip().upper()
    max_c = args.max_clusters
    receipt_lines = []

    with structure.StructureReader(str(in_path)) as reader:
        for i, st in enumerate(reader, 1):
            if max_c > 0 and i > max_c:
                break

            pose_path = "NO_POSE"
            if lig_name:
                lig_atoms = [
                    a.index for a in st.atom if a.pdbres.strip().upper() == lig_name
                ]
                if lig_atoms:
                    lig_st = st.extract(lig_atoms)
                    lig_st.title = lig_name
                    lig_st.property["s_m_title"] = lig_name
                    sdf_out = f"{args.prefix}_pose_{i}.sdf"
                    mae_out = f"{args.prefix}_pose_{i}.mae"
                    lig_st.write(sdf_out)
                    lig_st.write(mae_out)
                    st.deleteAtoms(lig_atoms)
                    pose_path = str(Path(sdf_out).resolve())

            c_title = f"{args.prefix}_centroid_{i}"
            st.title = c_title
            st.property["s_m_title"] = c_title
            pdb_out = f"{args.prefix}_centroid_{i}.pdb"
            mae_out = f"{args.prefix}_centroid_{i}.mae"
            st.write(pdb_out)
            st.write(mae_out)

            row = f"{c_title}\t{Path(pdb_out).resolve()}\t{Path(mae_out).resolve()}\t{pose_path}\n"
            receipt_lines.append(row)

    with open(args.out_receipt, "w", encoding="utf-8") as rf:
        rf.writelines(receipt_lines)

    print(f"Extracted {len(receipt_lines)} cluster centroid(s) -> {args.out_receipt}")


if __name__ == "__main__":
    main()
