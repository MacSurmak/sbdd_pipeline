"""
Common parsers for Glide CSV, Prime CSV, and GNINA SDF outputs.
"""

from pathlib import Path
import csv
import re
from typing import Dict, Any


def parse_glide_csv(csv_path: Path) -> Dict[str, Dict[str, float]]:
    data: Dict[str, Dict[str, float]] = {}
    counts: Dict[str, int] = {}
    if not csv_path.exists():
        return data
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw_t = (
                row.get("title")
                if row.get("title") is not None
                else (
                    row.get("Title")
                    if row.get("Title") is not None
                    else row.get("s_m_title")
                )
            )
            title = raw_t.strip() if raw_t else ""
            if not title:
                continue

            score: float | None = None
            for k in ["r_i_docking_score", "docking_score"]:
                if k in row and row[k]:
                    try:
                        score = float(row[k])
                        break
                    except ValueError as e:
                        raise ValueError(
                            f"FATAL: Unparseable GlideScore '{row[k]}' for ligand '{title}': {e}"
                        ) from e
            if score is None:
                for k in ["r_i_glide_gscore", "gscore"]:
                    if k in row and row[k]:
                        try:
                            score = float(row[k])
                            break
                        except ValueError as e:
                            raise ValueError(
                                f"FATAL: Unparseable GScore '{row[k]}' for ligand '{title}': {e}"
                            ) from e
            if score is not None:
                counts[title] = (counts[title] + 1) if title in counts else 1
                if title not in data or score < data[title]["GlideScore"]:
                    data[title] = {"GlideScore": score, "Glide_Modes": float(counts[title])}
                else:
                    data[title]["Glide_Modes"] = float(counts[title])
    return data


def parse_prime_csv(csv_path: Path) -> Dict[str, Dict[str, float]]:
    data: Dict[str, Dict[str, float]] = {}
    if not csv_path.exists():
        return data
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for row in reader:
            title = (row.get("title") or row.get("Title") or row.get("s_m_title") or "").strip()
            if not title:
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
                if title not in data or props.get("MMGBSA", float("inf")) < data[title].get("MMGBSA", float("inf")):
                    data[title] = props
    return data


def parse_gnina_sdf(sdf_path: Path, prefix: str = "") -> Dict[str, Dict[str, float]]:
    """
    Parses docked or GNINA-rescored SDF files.
    Selects the #1 representative pose per substance strictly by its native physical energy
    (GlideScore if present, otherwise Vina affinity) while attaching GNINA CNNscore and CNNaffinity
    as orthogonal metrics on that representative pose.
    """
    data: Dict[str, Dict[str, float]] = {}
    best_native_energy: Dict[str, float] = {}
    mode_counts: Dict[str, int] = {}
    if not sdf_path or not sdf_path.exists():
        return data

    from rdkit import Chem

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    for idx, mol in enumerate(suppl, start=1):
        if mol is None:
            raise ValueError(
                f"FATAL: Corrupted molecule record at index {idx} in {sdf_path}"
            )

        title = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        if not title and mol.HasProp("s_m_title"):
            title = mol.GetProp("s_m_title").strip()

        props: Dict[str, float] = {}
        native_glide: float | None = None
        native_vina: float | None = None

        for prop_name in mol.GetPropNames():
            tag_clean = prop_name.strip().lower()
            raw_val = mol.GetProp(prop_name).strip()
            if not raw_val or raw_val.lower() == "nan":
                continue

            if tag_clean in ("r_i_docking_score", "docking_score", "r_i_glide_gscore"):
                native_glide = float(raw_val)
            elif tag_clean == "unidock_vina_score":
                native_vina = float(raw_val.split()[0])
                props[f"{prefix}Vina"] = native_vina
            elif tag_clean in ("cnnscore", "cnnaffinity", "minimizedaffinity", "affinity"):
                try:
                    val = float(raw_val.split()[0])
                except ValueError as e:
                    raise ValueError(
                        f"FATAL: Unparseable metric '{prop_name}={raw_val}' for ligand '{title}' in {sdf_path}: {e}"
                    )
                if tag_clean == "cnnscore":
                    props[f"{prefix}CNNscore"] = val
                elif tag_clean == "cnnaffinity":
                    props[f"{prefix}CNNaffinity"] = val
                elif tag_clean in ("minimizedaffinity", "affinity"):
                    if f"{prefix}Vina" not in props:
                        props[f"{prefix}Vina"] = val
                    if native_vina is None:
                        native_vina = val

        if props and title:
            mode_counts[title] = (mode_counts[title] + 1) if title in mode_counts else 1
            primary_e = (
                native_glide
                if native_glide is not None
                else (native_vina if native_vina is not None else float("inf"))
            )
            prev_best = (
                best_native_energy[title]
                if title in best_native_energy
                else float("inf")
            )
            if title not in data or primary_e < prev_best:
                best_native_energy[title] = primary_e
                props[f"{prefix}Modes"] = float(mode_counts[title])
                data[title] = props
            else:
                data[title][f"{prefix}Modes"] = float(mode_counts[title])

    return data
