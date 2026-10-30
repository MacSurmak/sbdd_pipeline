#!/usr/bin/env python3
"""
bin/batch_dock_ligands.py
Batches prepared ligand SDF and MAEGZ files into chunks of size batch_size
for parallel Glide and GNINA execution on Slurm.
"""

import argparse
import gzip
import json
import shutil
import sys
from pathlib import Path


def merge_sdf_files(sdf_paths: list[Path], out_sdf: Path) -> None:
    with open(out_sdf, "w", encoding="utf-8") as out_f:
        for p in sdf_paths:
            if not p.exists() or p.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Prepared ligand SDF missing or empty: {p}"
                )
            content = p.read_text(encoding="utf-8", errors="ignore").strip()
            if not content:
                raise ValueError(f"FATAL: Prepared ligand SDF is empty: {p}")
            if not content.endswith("$$$$"):
                content += "\n$$$$"
            out_f.write(content + "\n")


def merge_maegz_files(mae_paths: list[Path], out_maegz: Path) -> None:
    with gzip.open(out_maegz, "wb") as out_gz:
        for p in mae_paths:
            if not p.exists() or p.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Prepared ligand MAEGZ missing or empty: {p}"
                )
            if p.suffix.lower() == ".maegz":
                with gzip.open(p, "rb") as in_gz:
                    shutil.copyfileobj(in_gz, out_gz)
            else:
                with open(p, "rb") as in_f:
                    shutil.copyfileobj(in_f, out_gz)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch prepared ligands into SDF and MAEGZ chunks for docking"
    )
    parser.add_argument("--dock_plan", required=True, help="Input dock_plan.json")
    parser.add_argument(
        "--batch_size", type=int, required=True, help="Number of substances per chunk"
    )
    parser.add_argument(
        "--out_manifest",
        default="chunks_manifest.json",
        help="Output JSON manifest of generated chunks",
    )
    args = parser.parse_args()

    if args.batch_size <= 0:
        print(
            f"FATAL: --batch_size must be a positive integer, got {args.batch_size}",
            file=sys.stderr,
        )
        sys.exit(1)

    plan_path = Path(args.dock_plan).resolve()
    if not plan_path.exists():
        print(f"FATAL: Dock plan not found: {plan_path}", file=sys.stderr)
        sys.exit(1)

    plan_data = json.loads(plan_path.read_text(encoding="utf-8"))
    mode = plan_data.get("mode", "matrix")
    ligands = plan_data.get("ligands", [])
    if not ligands:
        print("FATAL: No ligands specified in dock_plan.json!", file=sys.stderr)
        sys.exit(1)

    manifest = []

    if mode == "self":
        for idx, lig in enumerate(ligands):
            sdf_in = Path(lig["path_sdf"]).resolve()
            mae_in = Path(lig["path_mae"]).resolve()
            out_sdf = Path(f"chunk_{idx:03d}.sdf")
            out_mae = Path(f"chunk_{idx:03d}.maegz")
            merge_sdf_files([sdf_in], out_sdf)
            merge_maegz_files([mae_in], out_mae)
            manifest.append(
                {
                    "chunk_idx": idx,
                    "receptor": lig.get("receptor", ""),
                    "sdf": out_sdf.name,
                    "maegz": out_mae.name,
                    "substances": [lig.get("ligand_title", lig["name"])],
                }
            )
    else:
        b_size = int(args.batch_size)
        batches = [ligands[i : i + b_size] for i in range(0, len(ligands), b_size)]
        for idx, batch in enumerate(batches):
            sdf_list = [Path(item["path_sdf"]).resolve() for item in batch]
            mae_list = [Path(item["path_mae"]).resolve() for item in batch]
            out_sdf = Path(f"chunk_{idx:03d}.sdf")
            out_mae = Path(f"chunk_{idx:03d}.maegz")
            merge_sdf_files(sdf_list, out_sdf)
            merge_maegz_files(mae_list, out_mae)
            manifest.append(
                {
                    "chunk_idx": idx,
                    "receptor": "",
                    "sdf": out_sdf.name,
                    "maegz": out_mae.name,
                    "substances": [
                        item.get("ligand_title", item["name"]) for item in batch
                    ],
                }
            )

    Path(args.out_manifest).write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(
        f"Batched {len(ligands)} substance(s) into {len(manifest)} chunk(s) (batch_size={args.batch_size}, mode={mode})."
    )


if __name__ == "__main__":
    main()
