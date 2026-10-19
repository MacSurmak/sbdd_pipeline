#!/usr/bin/env python3
import argparse
from pathlib import Path


def cmd_merge_sdf(args):
    with open(args.out_sdf, "w") as out:
        for sdf in args.in_sdfs:
            lines = Path(sdf).read_text(errors="ignore").splitlines()
            # Strip trailing empty lines or stray $$$$ to prevent empty blocks
            while lines and (lines[-1].strip() == "" or lines[-1].strip() == "$$$$"):
                lines.pop()
            if not lines:
                continue

            # Sanitise RDKit dimensionality warnings by enforcing 3D tags
            for i in range(len(lines)):
                if i == 1 or (i > 1 and lines[i - 2] == "$$$$"):
                    if len(lines[i]) >= 22 and lines[i][20:22] == "2D":
                        lines[i] = lines[i][:20] + "3D" + lines[i][22:]
                    elif "2D" in lines[i]:
                        lines[i] = lines[i].replace("2D", "3D")

            out.write("\n".join(lines) + "\n$$$$\n")


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    p_merge = subparsers.add_parser("merge_sdf")
    p_merge.add_argument("--in_sdfs", nargs="+", required=True)
    p_merge.add_argument("--out_sdf", required=True)

    args = parser.parse_args()
    if args.subcommand == "merge_sdf":
        cmd_merge_sdf(args)


if __name__ == "__main__":
    main()
