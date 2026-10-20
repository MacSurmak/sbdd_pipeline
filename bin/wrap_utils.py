#!/usr/bin/env python3
import argparse
from pathlib import Path


def cmd_merge_sdf(args):
    from rdkit import Chem
    from rdkit import RDLogger

    # Suppress warnings during parsing of chunks (e.g. 2D tags from GNINA)
    RDLogger.DisableLog("rdApp.*")

    with Chem.SDWriter(args.out_sdf) as writer:
        for sdf in args.in_sdfs:
            suppl = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)
            for mol in suppl:
                if mol is not None:
                    writer.write(mol)


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
