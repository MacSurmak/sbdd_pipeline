#!/usr/bin/env python3
"""
bin/assemble_mae.py
Merges Apo receptor MAE with heterogen MAEs (cofactors/ions) to rebuild a Holo receptor for Glide.
Preserves Schrodinger-specific formal charges, partial charges, and hydrogen bond networks
established during PrepWizard.
"""
import argparse
from schrodinger import structure

def main():
    parser = argparse.ArgumentParser(description="Assemble Holo MAE for Glide")
    parser.add_argument("--rec_mae", required=True, help="Input pure Apo receptor MAE")
    parser.add_argument("--hetero_maes", nargs="*", default=[], help="Input heterogen MAEs (cofactors/ions)")
    parser.add_argument("--out_mae", required=True, help="Output assembled Holo MAE")
    args = parser.parse_args()

    st = structure.StructureReader.read(args.rec_mae)
    for het_path in args.hetero_maes:
        het_st = structure.StructureReader.read(het_path)
        st.extend(het_st)

    st.write(args.out_mae)

if __name__ == "__main__":
    main()
