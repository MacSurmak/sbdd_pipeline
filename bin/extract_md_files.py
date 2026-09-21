#!/usr/bin/env python3
"""
File: bin/extract_md_files.py
Description: Extracts receptor PDB and ligand MOL2 for MD setup.
Computes formal charge for Antechamber. Supports both Schrödinger API and RDKit fallback.
"""
import argparse
import sys
from pathlib import Path

def extract_via_schrodinger(lig_file, lig_title, rec_file, out_dir):
    from schrodinger import structure
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Receptor PDB
    rec_st = structure.StructureReader.read(rec_file)
    rec_out = out_dir / "receptor.pdb"
    structure.StructureWriter.write(rec_st, str(rec_out))

    # 2. Ligand
    charge = 0
    lig_found = False
    all_ligs = []
    
    with structure.StructureReader(lig_file) as reader:
        for i, st in enumerate(reader):
            if i == 0 and (st.property.get("b_glide_receptor", 0) == 1 or st.atom_total > 500):
                continue
            all_ligs.append(st)
            title = (st.property.get("s_m_title") or st.title or "").strip()
            if not lig_title or title == lig_title.strip() or title.lower() == lig_title.strip().lower():
                lig_mol2 = out_dir / "ligand.mol2"
                structure.StructureWriter.write(st, str(lig_mol2))
                charge = st.formal_charge
                lig_found = True
                break

    if not lig_found and len(all_ligs) == 1:
        st = all_ligs[0]
        lig_mol2 = out_dir / "ligand.mol2"
        structure.StructureWriter.write(st, str(lig_mol2))
        charge = st.formal_charge
        lig_found = True

    if not lig_found:
        raise ValueError(f"Ligand '{lig_title}' not found in {lig_file}")

    with open(out_dir / "lig_charge.txt", "w") as f:
        f.write(str(charge))
    print(f"Extracted files successfully. Ligand formal charge: {charge}")

def extract_via_rdkit(lig_file, lig_title, rec_file, out_dir):
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
        if mol is None: continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        if not lig_title or name == lig_title.strip() or name.lower() == lig_title.strip().lower():
            target_mol = mol
            break
            
    if target_mol is None and len(suppl) == 1 and suppl[0] is not None:
        target_mol = suppl[0]
        
    if target_mol is None:
        raise ValueError(f"Ligand '{lig_title}' not found in {lig_file}")
        
    charge = Chem.GetFormalCharge(target_mol)
    
    from rdkit.Chem import AllChem
    target_mol = Chem.AddHs(target_mol, addCoords=True)
    if not target_mol.GetNumConformers() or all(target_mol.GetConformer().GetAtomPosition(i).z == 0 for i in range(target_mol.GetNumAtoms())):
        AllChem.EmbedMolecule(target_mol)

    # Save mol2 via Chem.MolToMol2Block or openbabel / antechamber prep
    sdf_out = out_dir / "ligand.sdf"
    writer = Chem.SDWriter(str(sdf_out))
    writer.write(target_mol)
    writer.close()
    
    # Use antechamber to convert sdf to mol2
    import subprocess
    cmd = ["antechamber", "-i", str(sdf_out), "-fi", "sdf", "-o", str(out_dir / "ligand.mol2"), "-fo", "mol2"]
    import shutil
    antechamber_bin = shutil.which("antechamber") or str(Path(sys.prefix) / "bin" / "antechamber")
    cmd = [antechamber_bin, "-i", "ligand.sdf", "-fi", "sdf", "-o", "ligand.mol2", "-fo", "mol2"]
    subprocess.run(cmd, check=True, cwd=str(out_dir))
    # Clean up antechamber intermediate files
    for tmp_ac in out_dir.glob("ANTECHAMBER*"):
        tmp_ac.unlink(missing_ok=True)
    for tmp_info in out_dir.glob("ATOMTYPE*"):
        tmp_info.unlink(missing_ok=True)

    with open(out_dir / "lig_charge.txt", "w") as f:
        f.write(str(charge))
    print(f"Extracted files via RDKit. Formal charge: {charge}")

def main():
    parser = argparse.ArgumentParser(description="Extracts receptor and ligand files for MD")
    parser.add_argument("--lig_file", required=True, help="Input docked pose file (.maegz or .sdf)")
    parser.add_argument("--lig_title", default="", help="Title or ID of target ligand")
    parser.add_argument("--rec_file", required=True, help="Input receptor file (.mae or .pdb)")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    args = parser.parse_args()

    try:
        extract_via_schrodinger(args.lig_file, args.lig_title, args.rec_file, args.out_dir)
    except ImportError:
        extract_via_rdkit(args.lig_file, args.lig_title, args.rec_file, args.out_dir)
    except Exception as e:
        print(f"Error in extract_md_files: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

