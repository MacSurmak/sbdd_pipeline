"""
Common parsers for Glide CSV, Prime CSV, and GNINA SDF outputs.
"""

from pathlib import Path
import csv
import re
from typing import Dict, Any


def parse_glide_csv(csv_path: Path) -> Dict[str, Dict[str, float]]:
    data: Dict[str, Dict[str, float]] = {}
    if not csv_path.exists():
        return data
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("title") or row.get("Title") or row.get("s_m_title")
            if not title or title in data:
                continue

            score = None
            for k in ["r_i_docking_score", "docking_score"]:
                if k in row and row[k]:
                    try:
                        score = float(row[k])
                        break
                    except ValueError as e:
                        raise ValueError(
                            f"FATAL: Unparseable GlideScore '{row[k]}' for ligand '{title}': {e}"
                        )
            if score is None:
                for k in ["r_i_glide_gscore", "gscore"]:
                    if k in row and row[k]:
                        try:
                            score = float(row[k])
                            break
                        except ValueError as e:
                            raise ValueError(
                                f"FATAL: Unparseable GScore '{row[k]}' for ligand '{title}': {e}"
                            )
            if score is not None:
                data[title] = {"GlideScore": score}
    return data


def parse_prime_csv(csv_path: Path) -> Dict[str, Dict[str, float]]:
    data: Dict[str, Dict[str, float]] = {}
    if not csv_path.exists():
        return data
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = row.get("title") or row.get("Title") or row.get("s_m_title")
            if not title or title in data:
                continue
            props: Dict[str, float] = {}
            for k, v in row.items():
                if not k or not v:
                    continue
                k_clean = k.strip()
                if k_clean in ["r_psp_MMGBSA_dG_Bind", "MMGBSA dG Bind"]:
                    try:
                        props["MMGBSA"] = float(v)
                        break
                    except ValueError as e:
                        raise ValueError(
                            f"FATAL: Unparseable MMGBSA score '{v}' for ligand '{title}': {e}"
                        )

            for k, v in row.items():
                if not k or not v:
                    continue
                if k.strip() in ["r_psp_Lig_Strain_Energy", "Lig Strain Energy"]:
                    try:
                        props["Lig_Strain"] = float(v)
                        break
                    except ValueError as e:
                        raise ValueError(
                            f"FATAL: Unparseable Lig_Strain score '{v}' for ligand '{title}': {e}"
                        )

            if props:
                data[title] = props
    return data


def parse_gnina_sdf(sdf_path: Path, prefix: str = "") -> Dict[str, Dict[str, float]]:
    data: Dict[str, Dict[str, float]] = {}
    if not sdf_path or not sdf_path.exists():
        return data
    with open(sdf_path, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    for entry in content.split("$$$$"):
        if entry.startswith("\n"):
            entry = entry[1:]
        elif entry.startswith("\r\n"):
            entry = entry[2:]

        lines = entry.splitlines()
        if len(lines) < 4:
            continue

        raw_title = lines[0].strip()
        prop_title = None
        props: Dict[str, float] = {}

        for i, line in enumerate(lines):
            line_str = line.strip()
            if line_str.startswith("> <") and line_str.endswith(">"):
                tag_clean = line_str[3:-1].strip().lower()
                if i + 1 < len(lines):
                    val_str = lines[i + 1].strip()
                    if tag_clean in ["_name", "s_m_title", "title", "compound_id"]:
                        prop_title = val_str
                    elif tag_clean in [
                        "cnnscore",
                        "cnnaffinity",
                        "minimizedaffinity",
                        "affinity",
                    ]:
                        try:
                            val = float(val_str)
                            if tag_clean == "cnnscore":
                                props[f"{prefix}CNNscore"] = val
                            elif tag_clean == "cnnaffinity":
                                props[f"{prefix}CNNaffinity"] = val
                            elif tag_clean in ["minimizedaffinity", "affinity"]:
                                props[f"{prefix}VinaAff"] = val
                        except ValueError as e:
                            raise ValueError(
                                f"FATAL: Unparseable GNINA metric '{val_str}' for tag '{tag_clean}': {e}"
                            )

        title = prop_title or raw_title
        if not title:
            continue

        title_clean = re.sub(r"(_model_\d+|_dock_\d+|_\d+)$", "", title).strip()

        if props:
            if title not in data:
                data[title] = props
            if title_clean and title_clean not in data:
                data[title_clean] = props
    return data
