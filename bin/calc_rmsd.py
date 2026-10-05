"""
Worker script to calculate symmetry-aware heavy-atom RMSD between reference and docked poses.
"""

import argparse
import csv
import logging
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import rmsd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def calc_via_schrodinger(ref_path: Path, poses_path: Path, receptor_name: str, target_title: str) -> list[dict]:
    from schrodinger import structure
    from schrodinger.structutils import rmsd

    ref_st = structure.StructureReader.read(str(ref_path))
    results = []

    with structure.StructureReader(str(poses_path)) as reader:
        for i, st in enumerate(reader):
            if i == 0 and (st.property.get("b_glide_receptor", 0) == 1 or st.atom_total > 500):
                continue
            title = (st.property.get("s_m_title") or st.title or "").strip()
            if target_title and title != target_title:
                continue

            conf_eval = rmsd.ConformerRmsd(
                ref_st, st, asl_expr="NOT atom.element H", in_place=True
            )
            conf_eval.use_symmetry = True
            rmsd_val = round(float(conf_eval.calculate()), 3)

            results.append({
                "Receptor": receptor_name,
                "Title": title or target_title,
                "InSitu_RMSD_A": rmsd_val,
                "RMSD_Heavy": rmsd_val
            })
            break
    return results


def calc_via_rdkit(ref_path: Path, poses_path: Path, receptor_name: str, target_title: str) -> list[dict]:
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    def _read_mols(path: Path):
        mols = []
        if path.suffix.lower() == ".sdf":
            suppl = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
            for m in suppl:
                if m is not None:
                    mols.append(m)
        elif path.suffix.lower() == ".pdb":
            m = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
            if m:
                mols.append(m)
        return mols

    ref_mols = _read_mols(ref_path)
    if not ref_mols:
        raise ValueError(f"FATAL: Could not read reference structure from {ref_path}")
    m_ref = Chem.RemoveHs(ref_mols[0])

    docked_mols = _read_mols(poses_path)
    results = []

    for st in docked_mols:
        title = st.GetProp("_Name").strip() if st.HasProp("_Name") else ""
        if target_title and title and title != target_title:
            continue
        m_dock = Chem.RemoveHs(st)
        if m_dock.GetNumAtoms() != m_ref.GetNumAtoms():
            raise ValueError(
                f"FATAL: Heavy atom count mismatch between {ref_path.name} ({m_ref.GetNumAtoms()}) "
                f"and {poses_path.name} ({m_dock.GetNumAtoms()}) for compound '{title}'!"
            )
        val = round(float(rdMolAlign.GetBestRMS(m_dock, m_ref)), 3)
        results.append({
            "Receptor": receptor_name,
            "Title": title or target_title,
            "InSitu_RMSD_A": val,
            "RMSD_Heavy": val
        })
        break

    return results


def main():
    parser = argparse.ArgumentParser(description="Calculate In-Situ Heavy-Atom RMSD")
    parser.add_argument("--reference", required=True, help="Path to reference ligand structure (.sdf or .mae)")
    parser.add_argument("--poses", required=True, help="Path to docked ligand structure (.sdf or .maegz)")
    parser.add_argument("--receptor", default="", help="Receptor frame name")
    parser.add_argument("--title", default="", help="Ligand title")
    parser.add_argument("--output", required=True, help="Path to output CSV")
    args = parser.parse_args()

    ref_path = Path(args.reference).resolve()
    poses_path = Path(args.poses).resolve()
    out_path = Path(args.output).resolve()

    if not ref_path.exists():
        raise FileNotFoundError(f"FATAL: Reference ligand file {ref_path} does not exist!")
    if not poses_path.exists():
        raise FileNotFoundError(f"FATAL: Docked poses file {poses_path} does not exist!")

    results = []
    try:
        results = calc_via_schrodinger(ref_path, poses_path, args.receptor, args.title)
    except Exception as e:
        logging.info(f"Schrodinger RMSD unavailable ({e}), trying RDKit fallback...")
        results = calc_via_rdkit(ref_path, poses_path, args.receptor, args.title)

    if not results:
        raise ValueError(f"FATAL: RMSD calculation produced no results for {args.receptor} / {args.title}")

    fieldnames = ["Receptor", "Title", "InSitu_RMSD_A", "RMSD_Heavy"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    logging.info(f"Successfully calculated in-situ RMSD: {results[0]['InSitu_RMSD_A']} Å -> {out_path}")


if __name__ == "__main__":
    main()
