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
        description="Batch prepared ligands into independent SDF (GNINA) and MAEGZ (Glide) chunks"
    )
    parser.add_argument("--dock_plan", required=True, help="Input dock_plan.json")
    parser.add_argument(
        "--glide_batch_size",
        type=int,
        required=True,
        help="Substances per Glide chunk (0 if Glide disabled)",
    )
    parser.add_argument(
        "--gnina_batch_size",
        type=int,
        required=True,
        help="Substances per GNINA chunk (0 if GNINA disabled)",
    )
    parser.add_argument(
        "--out_manifest",
        default="chunks_manifest.json",
        help="Output JSON manifest of generated chunks",
    )
    args = parser.parse_args()

    if args.glide_batch_size <= 0 and args.gnina_batch_size <= 0:
        print(
            "FATAL: At least one of --glide_batch_size or --gnina_batch_size must be positive!",
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

    glide_chunks = []
    gnina_chunks = []

    if mode == "self":
        for idx, lig in enumerate(ligands):
            if args.glide_batch_size > 0:
                mae_in = Path(lig["path_mae"]).resolve()
                out_mae = Path(f"glide_chunk_{idx:03d}.maegz")
                merge_maegz_files([mae_in], out_mae)
                glide_chunks.append(
                    {
                        "chunk_idx": idx,
                        "receptor": lig.get("receptor", ""),
                        "file": out_mae.name,
                        "substances": [lig.get("ligand_title", lig["name"])],
                    }
                )
            if args.gnina_batch_size > 0:
                sdf_in = Path(lig["path_sdf"]).resolve()
                out_sdf = Path(f"gnina_chunk_{idx:03d}.sdf")
                merge_sdf_files([sdf_in], out_sdf)
                gnina_chunks.append(
                    {
                        "chunk_idx": idx,
                        "receptor": lig.get("receptor", ""),
                        "file": out_sdf.name,
                        "substances": [lig.get("ligand_title", lig["name"])],
                    }
                )
    else:
        if args.glide_batch_size > 0:
            g_size = int(args.glide_batch_size)
            num_g_chunks = max(1, (len(ligands) + g_size - 1) // g_size)
            g_batches = [ligands[i::num_g_chunks] for i in range(num_g_chunks)]
            for idx, batch in enumerate(g_batches):
                mae_list = [Path(item["path_mae"]).resolve() for item in batch]
                out_mae = Path(f"glide_chunk_{idx:03d}.maegz")
                merge_maegz_files(mae_list, out_mae)
                glide_chunks.append(
                    {
                        "chunk_idx": idx,
                        "receptor": "",
                        "file": out_mae.name,
                        "substances": [
                            item.get("ligand_title", item["name"]) for item in batch
                        ],
                    }
                )

        if args.gnina_batch_size > 0:
            n_size = int(args.gnina_batch_size)
            num_n_chunks = max(1, (len(ligands) + n_size - 1) // n_size)
            n_batches = [ligands[i::num_n_chunks] for i in range(num_n_chunks)]
            for idx, batch in enumerate(n_batches):
                sdf_list = [Path(item["path_sdf"]).resolve() for item in batch]
                out_sdf = Path(f"gnina_chunk_{idx:03d}.sdf")
                merge_sdf_files(sdf_list, out_sdf)
                gnina_chunks.append(
                    {
                        "chunk_idx": idx,
                        "receptor": "",
                        "file": out_sdf.name,
                        "substances": [
                            item.get("ligand_title", item["name"]) for item in batch
                        ],
                    }
                )

    manifest = {"glide_chunks": glide_chunks, "gnina_chunks": gnina_chunks}
    Path(args.out_manifest).write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(
        f"Batched {len(ligands)} substance(s): {len(glide_chunks)} Glide chunk(s) (size={args.glide_batch_size}), "
        f"{len(gnina_chunks)} GNINA chunk(s) (size={args.gnina_batch_size}) [mode={mode}]."
    )


if __name__ == "__main__":
    main()
