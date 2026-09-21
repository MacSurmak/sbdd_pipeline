#!/usr/bin/env python3
"""
bin/split_sdf.py
Splits an SDF library into chunks of N molecules.
"""
import argparse
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Split SDF into batches")
    parser.add_argument("--input", required=True, help="Input SDF file")
    parser.add_argument("--batch_size", "--chunk_size", dest="batch_size", type=int, default=15, help="Molecules per chunk")
    parser.add_argument("--outdir", default=".", help="Output directory for chunks")
    args = parser.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.input, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    mols = [m.strip() for m in content.split("$$$$") if m.strip()]
    print(f"Detected {len(mols)} molecules in {args.input}")

    total_chunks = 0
    for idx, i in enumerate(range(0, len(mols), args.batch_size)):
        chunk = mols[i:i + args.batch_size]
        out_file = out_dir / f"chunk_{idx:03d}.sdf"
        with open(out_file, "w", encoding="utf-8") as out:
            out.write("\n$$$$\n".join(chunk) + "\n$$$$\n")
        total_chunks += 1

    print(f"Successfully created {total_chunks} chunks of size <= {args.batch_size}.")

if __name__ == "__main__":
    main()
