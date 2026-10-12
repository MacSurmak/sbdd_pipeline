#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--healed", required=True)
    parser.add_argument("--original", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    healed_lines = []
    with open(args.healed, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "TER")):
                healed_lines.append(line)
    
    het_lines = []
    with open(args.original, 'r') as f:
        for line in f:
            if line.startswith("HETATM"):
                # Exclude water if needed, but safe to keep all original HETATM 
                # (water clashing will be handled by PrepWizard/Declash downstream)
                het_lines.append(line)

    with open(args.out, 'w') as f:
        f.writelines(healed_lines)
        f.writelines(het_lines)
        f.write("END\n")

if __name__ == "__main__":
    main()
