"""
Worker script to calculate symmetry-aware heavy-atom RMSD between reference and docked poses.
"""
import argparse
import csv
import logging
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import rmsd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def main():
    parser = argparse.ArgumentParser(description="Calculate Conformer RMSD")
    parser.add_argument("--reference", required=True, help="Path to reference ligand MAE")
    parser.add_argument("--poses", required=True, help="Path to Glide PV or Ligand MAE")
    parser.add_argument("--output", required=True, help="Path to output CSV")
    args = parser.parse_args()

    ref_st = structure.StructureReader.read(args.reference)
    
    results = []
    
    with structure.StructureReader(args.poses) as reader:
        for i, st in enumerate(reader):
            # Skip the receptor
            if i == 0 and st.atom_total > 500:
                continue
                
            try:
                conf_eval = rmsd.ConformerRmsd(ref_st, st, asl_expr='NOT atom.element H', in_place=True)
                conf_eval.use_symmetry = True
                rmsd_val = conf_eval.calculate()
                
                gscore = st.property.get('r_i_glide_gscore', '')
                dscore = st.property.get('r_i_glide_docking_score', '')
                title = st.title
                
                # Extract GNINA energy and CNN properties if present in SDF
                gnina_affinity, cnn_score, cnn_affinity = "", "", ""
                for k, v in st.property.items():
                    k_low = k.lower()
                    if "minimizedaffinity" in k_low or (k_low.endswith("affinity") and "cnn" not in k_low):
                        gnina_affinity = v
                    elif "cnnscore" in k_low:
                        cnn_score = v
                    elif "cnnaffinity" in k_low:
                        cnn_affinity = v
                
                # Universal primary score (GlideScore or GNINA affinity)
                primary_score = gscore or gnina_affinity or dscore or ""

                results.append({
                    "Title": title,
                    "Pose_Index": i,
                    "Score": primary_score,
                    "GlideScore": gscore,
                    "DockingScore": dscore,
                    "CNN_Score": cnn_score,
                    "CNN_Affinity": cnn_affinity,
                    "RMSD_Heavy": round(rmsd_val, 3)
                })
            except Exception as e:
                logging.error(f"Error processing pose {i}: {e}")

    keys = ["Title", "Pose_Index", "Score", "GlideScore", "DockingScore", "CNN_Score", "CNN_Affinity", "RMSD_Heavy"]

    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)
        
    logging.info(f"Successfully calculated RMSD for {len(results)} poses. Saved to {args.output}")

if __name__ == "__main__":
    main()