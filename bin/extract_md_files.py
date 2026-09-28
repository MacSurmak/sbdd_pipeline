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

    # 2. Extract Ligand if provided
    if (
        args.lig_file
        and args.lig_file.lower() != "none"
        and Path(args.lig_file).exists()
    ):
        with structure.StructureReader(args.lig_file) as reader:
            for i, st in enumerate(reader):
                st.title = "UNK"
                for res in st.residue:
                    res.pdbres = "UNK"
                lig_out = org_dir / "UNK.mol2"
                st.write(str(lig_out))
                (org_dir / "UNK_charge.txt").write_text(str(st.formal_charge))
                organics.append(
                    {"name": "UNK", "mol2": "UNK.mol2", "charge": "UNK_charge.txt"}
                )
                break

    with open(out_dir / "organics.json", "w") as f:
        json.dump(organics, f)

    print(f"Extracted dry receptor and {len(organics)} organic molecule(s).")


def extract_via_rdkit(args):
    from rdkit import Chem

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Receptor copy / conversion
    rec_path = Path(rec_file)
    if rec_path.suffix.lower() == ".pdb":
        import shutil

        shutil.copy(rec_path, out_dir / "receptor.pdb")

    suppl = Chem.SDMolSupplier(str(lig_file), removeHs=False)
    target_mol = None

    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        if (
            not lig_title
            or name == lig_title.strip()
            or name.lower() == lig_title.strip().lower()
        ):
            target_mol = mol
            break

    if target_mol is None:
        raise ValueError(f"FATAL: Ligand '{lig_title}' not found in {lig_file}")

    charge = Chem.GetFormalCharge(target_mol)

    from rdkit.Chem import AllChem

    target_mol = Chem.AddHs(target_mol, addCoords=True)
    if not target_mol.GetNumConformers() or all(
        target_mol.GetConformer().GetAtomPosition(i).z == 0
        for i in range(target_mol.GetNumAtoms())
    ):
        AllChem.EmbedMolecule(target_mol)

    # Save mol2 via Chem.MolToMol2Block or openbabel / antechamber prep
    sdf_out = out_dir / "ligand.sdf"
    writer = Chem.SDWriter(str(sdf_out))
    writer.write(target_mol)
    writer.close()

    # Use antechamber to convert sdf to mol2 inside out_dir only
    import subprocess
    import shutil

    antechamber_bin = shutil.which("antechamber") or str(
        Path(sys.prefix) / "bin" / "antechamber"
    )
    cmd = [
        antechamber_bin,
        "-i",
        "ligand.sdf",
        "-fi",
        "sdf",
        "-o",
        "ligand.mol2",
        "-fo",
        "mol2",
    ]
    subprocess.run(cmd, check=True, cwd=str(out_dir))
    # Clean up antechamber intermediate files in out_dir
    for tmp_ac in out_dir.glob("ANTECHAMBER*"):
        tmp_ac.unlink(missing_ok=True)
    for tmp_info in out_dir.glob("ATOMTYPE*"):
        tmp_info.unlink(missing_ok=True)

    with open(out_dir / "lig_charge.txt", "w") as f:
        f.write(str(charge))
    print(f"Extracted files via RDKit. Formal charge: {charge}")


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
