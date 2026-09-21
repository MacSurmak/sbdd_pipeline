"""
Worker script executed via $SCHRODINGER/run python3.
Filters a Pose Viewer file, keeping only ligands explicitly present in the whitelist for a given frame.
"""
import argparse
import sys
import csv
from schrodinger import structure

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input _pv.maegz file")
    parser.add_argument("--output", required=True, help="Filtered output _pv.maegz file")
    parser.add_argument("--whitelist", required=True, help="Path to whitelist.csv")
    parser.add_argument("--frame", required=True, help="Current frame ID")
    parser.add_argument("--max_poses", type=int, default=1, help="Maximum poses per ligand for MM-GBSA")
    args = parser.parse_args()

    allowed_titles = set()
    with open(args.whitelist, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if str(row["Frame"]) == str(args.frame):
                allowed_titles.add(row["Title"])

    kept = 0
    pose_counts = {}
    with structure.StructureReader(args.input) as reader, structure.StructureWriter(args.output) as writer:
        for i, st in enumerate(reader):
            # Always keep the receptor (first structure)
            if i == 0 and (st.property.get("b_glide_receptor", 0) == 1 or st.atom_total > 500):
                writer.append(st)
                continue
                
            title = (st.property.get("s_m_title") or st.title or "").strip()
            # Only keep top n poses for whitelisted ligands
            if title in allowed_titles and pose_counts.get(title, 0) < args.max_poses:
                writer.append(st)
                pose_counts[title] = pose_counts.get(title, 0) + 1
                kept += 1

    print(f"SUCCESS: Extracted receptor and {kept} whitelisted poses for Frame {args.frame}.")

if __name__ == "__main__":
    main()