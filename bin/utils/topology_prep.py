"""
File: src/utils/topology_prep.py
Description: Cleaned and refactored logic for PDB merging and tLEaP bond detection.
Replaces the legacy 'combine_to_complex.py' and 'detect_bonds.py'.
"""
import math
from pathlib import Path

def merge_protein_ligand(prot_pdb: Path, lig_mol2: Path, out_pdb: Path) -> int:
    """
    Merges AMBER-prepped protein PDB with ligand Mol2.
    Returns the number of ligand atoms inserted.
    """
    lines = []
    last_serial = 0
    used_chains = set()

    with open(prot_pdb, "r") as f:
        for line in f:
            rec = line[0:6]
            if rec.startswith(("CONECT", "MASTER", "END")):
                continue
            if rec.startswith(("ATOM", "HETATM")):
                try:
                    s = int(line[6:11])
                    if s > last_serial: last_serial = s
                except ValueError: pass
                if len(line) >= 22 and line[21].strip():
                    used_chains.add(line[21])
            lines.append(line)

    ligand_chain = next((c for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ" if c not in used_chains), "Z")
    lig_atoms = []

    with open(lig_mol2, "r") as f:
        in_atom = False
        for line in f:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and not line.startswith("@<TRIPOS>ATOM"):
                in_atom = False
                continue
            if in_atom and line:
                parts = line.split()
                if len(parts) >= 6:
                    elem = parts[5].split(".")[0][:2].capitalize() if "." in parts[5] else parts[5][0].upper()
                    lig_atoms.append({
                        "name": parts[1][:4], "x": float(parts[2]), "y": float(parts[3]), 
                        "z": float(parts[4]), "elem": elem
                    })

    serial = last_serial
    ligand_lines = []
    for a in lig_atoms:
        serial += 1
        line = (
            f"HETATM{serial:>5} {a['name']:>4} UNK {ligand_chain}   1       "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00  0.00          {a['elem']:>2}  \n"
        )
        ligand_lines.append(line)

    with open(out_pdb, "w") as out:
        out.writelines(lines)
        out.write("TER\n")
        out.writelines(ligand_lines)
        out.write("END\n")

    return len(lig_atoms)

def generate_tleap_bonds(pdb_path: Path, unit_name: str = "complex") -> str:
    """
    Parses PDB for Disulfide bridges (SG-SG distance) and generates tLEaP bond commands.
    """
    cys_atoms = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                if line[12:16].strip() == "SG" and line[17:20].strip() in ["CYS", "CYX", "CYM"]:
                    cys_atoms.append({
                        "chain": line[21], "resSeq": int(line[22:26]),
                        "name": "SG", "x": float(line[30:38]),
                        "y": float(line[38:46]), "z": float(line[46:54])
                    })

    bonds = []
    n = len(cys_atoms)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = cys_atoms[i], cys_atoms[j]
            if a["resSeq"] == b["resSeq"] and a["chain"] == b["chain"]: continue
            dist = math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)
            if dist <= 2.3:
                bonds.append(f"bond {unit_name}.{a['resSeq']}.SG {unit_name}.{b['resSeq']}.SG")

    return "\n".join(bonds)