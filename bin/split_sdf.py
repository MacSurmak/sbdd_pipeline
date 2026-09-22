#!/usr/bin/env python3
"""
bin/split_sdf.py
Splits an SDF library into chunks of N molecules.
Splits an SDF library into chunks of N molecules with minimal memory footprint (streaming).
"""
import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Split SDF into batches")
    parser.add_argument("--input", required=True, help="Input SDF file")
    parser.add_argument(
        "--batch_size",
        "--chunk_size",
        dest="batch_size",
        type=int,
        default=15,
        help="Molecules per chunk",
    )
    parser.add_argument("--outdir", default=".", help="Output directory for chunks")
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    total_mols = 0
    chunk_idx = 0
    current_chunk = []
    current_mol_lines = []

    with open(args.input, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            current_mol_lines.append(line)
            if line.startswith("$$$$"):
                mol_text = "".join(current_mol_lines).strip()
                if mol_text:
                    current_chunk.append(mol_text)
                    total_mols += 1
                current_mol_lines = []
                if len(current_chunk) >= args.batch_size:
                    out_file = out_dir / f"chunk_{chunk_idx:03d}.sdf"
                    with open(out_file, "w", encoding="utf-8") as out:
                        out.write("\n$$$$\n".join(current_chunk) + "\n$$$$\n")
                    chunk_idx += 1
                    current_chunk = []

        if current_mol_lines:
            mol_text = "".join(current_mol_lines).strip()
            if mol_text:
                current_chunk.append(mol_text)
                total_mols += 1

        if current_chunk:
            out_file = out_dir / f"chunk_{chunk_idx:03d}.sdf"
            with open(out_file, "w", encoding="utf-8") as out:
                out.write("\n$$$$\n".join(current_chunk) + "\n$$$$\n")
            chunk_idx += 1

    print(
        f"Processed {total_mols} molecules across {chunk_idx} chunks (size <= {args.batch_size})."
    )


if __name__ == "__main__":
    main()
