#!/usr/bin/env python3
"""
File: bin/extract_md_files.py
Description: Extracts receptor PDB and ligand MOL2 for MD setup.
Computes formal charge for Antechamber. Supports both Schrödinger API and RDKit fallback.
"""
import argparse
import sys
from pathlib import Path


ION_RESNAMES = {"ZN", "MG", "CA", "NA", "K", "CL", "ZN2", "MG2", "CA2"}


def extract_via_schrodinger(args):
    from schrodinger import structure
    from schrodinger.structutils import analyze
    import json

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    org_dir = out_dir / "org_files"
    org_dir.mkdir(parents=True, exist_ok=True)

    rec_st = structure.StructureReader.read(args.rec_file)
    organics = []

    # 1. Extract organic cofactors (Leave ions in protein)
    if args.cofactors:
        for cof in args.cofactors:
            cof = cof.strip()
            if not cof or cof.upper() in ION_RESNAMES:
                continue
            asl = f"res.ptype '{cof}'"
            atoms = analyze.evaluate_asl(rec_st, asl)
            if atoms:
                cof_st = rec_st.extract(atoms)
                cof_st.title = cof
                for res in cof_st.residue:
                    res.pdbres = cof
                cof_out = org_dir / f"{cof}.mol2"
                cof_st.write(str(cof_out))
                (org_dir / f"{cof}_charge.txt").write_text(str(cof_st.formal_charge))
                organics.append(
                    {"name": cof, "mol2": f"{cof}.mol2", "charge": f"{cof}_charge.txt"}
                )
                rec_st.deleteAtoms(atoms)

    # Save dry receptor (Protein + Ions)
    rec_out = out_dir / "receptor_dry.pdb"
    rec_st.write(str(rec_out))

    # 2. Extract Ligand preserving canonical residue name
    lig_code = (args.lig_name or "LIG").strip()[:4].upper()
    if (
        args.lig_file
        and args.lig_file.lower() != "none"
        and Path(args.lig_file).exists()
    ):
        with structure.StructureReader(args.lig_file) as reader:
            for i, st in enumerate(reader):
                st.title = lig_code
                for res in st.residue:
                    res.pdbres = lig_code
                lig_out = org_dir / f"{lig_code}.mol2"
                st.write(str(lig_out))
                (org_dir / f"{lig_code}_charge.txt").write_text(str(st.formal_charge))
                organics.append(
                    {"name": lig_code, "mol2": f"{lig_code}.mol2", "charge": f"{lig_code}_charge.txt"}
                )
                break

    with open(out_dir / "organics.json", "w") as f:
        json.dump(organics, f)

    print(f"Extracted dry receptor and {len(organics)} organic molecule(s).")


def extract_via_rdkit(args):
    from rdkit import Chem
    import json

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    org_dir = out_dir / "org_files"
    org_dir.mkdir(parents=True, exist_ok=True)

    # 1. Receptor copy (simplistic, RDKit doesn't easily split cofactors without heavy processing)
    rec_path = Path(args.rec_file)
    if rec_path.suffix.lower() == ".pdb":
        import shutil

        shutil.copy(rec_path, out_dir / "receptor_dry.pdb")

    organics = []

    # 2. Extract Ligand if provided
    if (
        args.lig_file
        and args.lig_file.lower() != "none"
        and Path(args.lig_file).exists()
    ):
        suppl = Chem.SDMolSupplier(str(args.lig_file), removeHs=False)
        target_mol = next((m for m in suppl if m is not None), None)

        if target_mol is None:
            raise ValueError(f"FATAL: No valid ligand found in {args.lig_file}")

        charge = Chem.GetFormalCharge(target_mol)

        lig_code = (args.lig_name or "LIG").strip()[:4].upper()
        sdf_out = org_dir / f"{lig_code}.sdf"
        writer = Chem.SDWriter(str(sdf_out))
        writer.write(target_mol)
        writer.close()

        import subprocess
        import shutil

        antechamber_bin = shutil.which("antechamber") or str(
            Path(sys.prefix) / "bin" / "antechamber"
        )
        cmd = [
            antechamber_bin,
            "-i",
            f"{lig_code}.sdf",
            "-fi",
            "sdf",
            "-o",
            f"{lig_code}.mol2",
            "-fo",
            "mol2",
            "-rn",
            lig_code,
        ]
        subprocess.run(cmd, check=True, cwd=str(org_dir))

        for tmp in org_dir.glob("ANTECHAMBER*"):
            tmp.unlink(missing_ok=True)
        for tmp in org_dir.glob("ATOMTYPE*"):
            tmp.unlink(missing_ok=True)

        (org_dir / f"{lig_code}_charge.txt").write_text(str(charge))
        organics.append({"name": lig_code, "mol2": f"{lig_code}.mol2", "charge": f"{lig_code}_charge.txt"})

    with open(out_dir / "organics.json", "w") as f:
        json.dump(organics, f)

    print(f"Extracted dry receptor and {len(organics)} organic molecule(s) via RDKit.")


def main():
    parser = argparse.ArgumentParser(
        description="Extracts receptor and ligand files for MD"
    )
    parser.add_argument(
        "--rec_file", required=True, help="Input prepared receptor file (.mae or .pdb)"
    )
    parser.add_argument(
        "--lig_file", default="none", help="Input docked pose file (.sdf or .mae)"
    )
    parser.add_argument(
        "--lig_name", default="", help="Canonical 3-letter residue name for ligand"
    )
    parser.add_argument("--cofactors", nargs="*", default=[], help="Cofactors list")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    args = parser.parse_args()

    try:
        extract_via_schrodinger(args)
    except ImportError:
        extract_via_rdkit(args)
    except Exception as e:
        print(f"Error in extract_md_files: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
