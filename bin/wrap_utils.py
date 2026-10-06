#!/usr/bin/env python3
import argparse
from pathlib import Path


def cmd_merge_sdf(args):
    with open(args.out_sdf, "w") as out:
        for sdf in args.in_sdfs:
            content = Path(sdf).read_text(errors="ignore")
            out.write(content)
            if not content.endswith("$$$$\n"):
                out.write("\n$$$$\n")


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
