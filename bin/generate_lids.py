#!/usr/bin/env python3
"""
File: src/schrodinger_workers/generate_lids.py
Description: Worker script executed via $SCHRODINGER/run python3.
Generates publication-quality 2D Ligand Interaction Diagrams (LID) using
the Schrödinger LiveDesign API (schrodinger.application.livedesign.lid).
Enforces headless offscreen Qt initialization to prevent segmentation faults.
"""
import argparse
import json
import os
import re
import sys
from pathlib import Path

# 1. Enforce headless offscreen platform before any Qt or Schrodinger module imports
os.environ["QT_QPA_PLATFORM"] = "offscreen"

# 2. Initialize QApplication instance to satisfy C++ Qt paint engine in _canvas2d
try:
    from PyQt6.QtWidgets import QApplication
except ImportError:
    from PyQt5.QtWidgets import QApplication

app = QApplication.instance()
if app is None:
    app = QApplication(["--platform", "offscreen"])

from schrodinger import structure
from schrodinger.application.livedesign import lid


def _clean_title(title: str) -> str:
    """Strips conformer and docking mode suffixes for title matching."""
    return re.sub(r"(_model_\d+|_dock_\d+|_\d+)$", "", title).strip()


def _render_and_save(lig_st, task: dict, receptor_cache: dict) -> bool:
    """Renders 2D interaction diagram for a ligand-receptor pair and saves PNG."""
    rec_file = task.get("receptor_file")
    out_png = Path(task["out_png"])
    title = task["title"]

    rec_st = None
    if rec_file and Path(rec_file).exists():
        if rec_file not in receptor_cache:
            try:
                receptor_cache[rec_file] = structure.StructureReader.read(rec_file)
            except Exception as e:
                print(f"[ERROR] Failed reading receptor file {rec_file}: {e}", file=sys.stderr)
                receptor_cache[rec_file] = None
        rec_st = receptor_cache[rec_file]

    if rec_st is None:
        print(f"[ERROR] No receptor structure available for {title} (path: {rec_file})", file=sys.stderr)
        return False

    try:
        qimage = lid.generate_lid(lig_st=lig_st, rec_st=rec_st)
        if qimage is not None and not qimage.isNull():
            out_png.parent.mkdir(parents=True, exist_ok=True)
            qimage.save(str(out_png), "PNG")
            print(f"[SUCCESS] Saved LID: {out_png.name}")
            return True
        else:
            print(f"[WARN] Empty LID generated for compound: {title}", file=sys.stderr)
            return False
    except Exception as e:
        print(f"[ERROR] Rendering failed for {title}: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description="Schrodinger 2D LID Generator Worker")
    parser.add_argument("--manifest", required=True, help="Path to JSON manifest file")
    parser.add_argument("--out_dir", required=True, help="Directory to save generated PNGs")
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        print(f"[CRITICAL] Manifest file not found: {manifest_path}", file=sys.stderr)
        sys.exit(1)

    with open(manifest_path, "r", encoding="utf-8") as f:
        tasks = json.load(f)

    print(f"Loaded {len(tasks)} LID task(s) from manifest.")

    receptor_cache = {}
    tasks_by_file = {}
    for t in tasks:
        tasks_by_file.setdefault(t["ligand_file"], []).append(t)

    success_count = 0

    for lig_file_str, file_tasks in tasks_by_file.items():
        lig_path = Path(lig_file_str)
        if not lig_path.exists():
            print(f"[ERROR] Ligand source file missing: {lig_path}", file=sys.stderr)
            continue

        # Build lookup table for fast matching
        task_lookup = {}
        for t in file_tasks:
            raw_t = t["title"].strip()
            task_lookup[raw_t] = t
            task_lookup[_clean_title(raw_t)] = t
            task_lookup[raw_t.lower()] = t

        rendered_outputs = set()
        collected_ligands = []

        try:
            with structure.StructureReader(str(lig_path)) as reader:
                for idx, st in enumerate(reader):
                    # Skip receptor entry in Pose Viewer files
                    if idx == 0 and (st.property.get("b_glide_receptor", 0) == 1 or st.atom_total > 500):
                        continue

                    collected_ligands.append(st)
                    st_title = (
                        st.property.get("s_m_title")
                        or st.title
                        or st.property.get("_Name")
                        or f"pose_{idx}"
                    ).strip()

                    match = (
                        task_lookup.get(st_title)
                        or task_lookup.get(_clean_title(st_title))
                        or task_lookup.get(st_title.lower())
                    )

                    if match:
                        out_key = match["out_png"]
                        if out_key not in rendered_outputs:
                            if _render_and_save(st, match, receptor_cache):
                                success_count += 1
                            rendered_outputs.add(out_key)

            # Fallback for dedicated single-ligand files (e.g. ref_native_f*.mae)
            if len(file_tasks) == 1 and len(rendered_outputs) == 0 and collected_ligands:
                single_task = file_tasks[0]
                if _render_and_save(collected_ligands[0], single_task, receptor_cache):
                    success_count += 1

        except Exception as e:
            print(f"[ERROR] Failed reading ligand poses from {lig_path}: {e}", file=sys.stderr)

    print(f"LID generation finished. Successfully generated: {success_count}/{len(tasks)} diagrams.")


if __name__ == "__main__":
    main()