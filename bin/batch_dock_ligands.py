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
    with open(out_maegz, "wb") as out_gz:
        for p in mae_paths:
            if not p.exists() or p.stat().st_size == 0:
                raise FileNotFoundError(
                    f"FATAL: Prepared ligand MAEGZ missing or empty: {p}"
                )
            with open(p, "rb") as in_f:
                shutil.copyfileobj(in_f, out_gz)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch prepared ligands into independent SDF (Uni-Dock) and MAEGZ (Glide) chunks"
    )
    parser.add_argument("--dock_plan", required=True, help="Input dock_plan.json")
    parser.add_argument(
        "--glide_batch_size",
        type=int,
        required=True,
        help="Substances per Glide chunk (0 if Glide disabled)",
    )
    parser.add_argument(
        "--unidock_batch_size",
        type=int,
        required=True,
        help="Substances per Uni-Dock chunk (0 if Uni-Dock disabled)",
    )
    parser.add_argument(
        "--out_manifest",
        default="chunks_manifest.json",
        help="Output JSON manifest of generated chunks",
    )
    args = parser.parse_args()

    if args.glide_batch_size <= 0 and args.unidock_batch_size <= 0:
        print(
            "FATAL: At least one of --glide_batch_size or --unidock_batch_size must be positive!",
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
    unidock_chunks = []

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
            if args.unidock_batch_size > 0:
                sdf_in = Path(lig["path_sdf"]).resolve()
                out_sdf = Path(f"unidock_chunk_{idx:03d}.sdf")
                merge_sdf_files([sdf_in], out_sdf)
                unidock_chunks.append(
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

        if args.unidock_batch_size > 0:
            u_size = int(args.unidock_batch_size)
            num_u_chunks = max(1, (len(ligands) + u_size - 1) // u_size)
            u_batches = [ligands[i::num_u_chunks] for i in range(num_u_chunks)]
            for idx, batch in enumerate(u_batches):
                sdf_list = [Path(item["path_sdf"]).resolve() for item in batch]
                out_sdf = Path(f"unidock_chunk_{idx:03d}.sdf")
                merge_sdf_files(sdf_list, out_sdf)
                unidock_chunks.append(
                    {
                        "chunk_idx": idx,
                        "receptor": "",
                        "file": out_sdf.name,
                        "substances": [
                            item.get("ligand_title", item["name"]) for item in batch
                        ],
                    }
                )

    manifest = {"glide_chunks": glide_chunks, "unidock_chunks": unidock_chunks}
    Path(args.out_manifest).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(
        f"Batched {len(ligands)} substance(s): {len(glide_chunks)} Glide chunk(s) (size={args.glide_batch_size}), "
        f"{len(unidock_chunks)} Uni-Dock chunk(s) (size={args.unidock_batch_size}) [mode={mode}]."
    )


if __name__ == "__main__":
    main()
