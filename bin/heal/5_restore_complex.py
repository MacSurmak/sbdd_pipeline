#!/usr/bin/env python3
import argparse
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--healed", required=True)
    parser.add_argument("--hetatm", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    healed_lines = []
    with open(args.healed, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "TER")):
                healed_lines.append(line)

    het_lines = []
    with open(args.hetatm, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                het_lines.append(line)

    with open(args.out, "w", encoding="utf-8") as f:
        f.writelines(healed_lines)
        if healed_lines and not healed_lines[-1].startswith("TER"):
            f.write("TER\n")
        f.writelines(het_lines)
        f.write("END\n")


if __name__ == "__main__":
    main()
