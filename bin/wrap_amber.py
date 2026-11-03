#!/usr/bin/env python3
import argparse
import sys
import subprocess
from pathlib import Path


def run_cmd(cmd):
    print(f"Exec: {' '.join(cmd)}")
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        print(f"FATAL: Command failed: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(res.returncode)


def parse_mol2_resname(mol2_file: Path) -> str:
    """Extracts strict 3-4 letter residue name from canonical Tripos substructure record."""
    in_sub = False
    for line in mol2_file.read_text(errors="ignore").splitlines():
        if line.startswith("@<TRIPOS>SUBSTRUCTURE"):
            in_sub = True
            continue
        if in_sub and line.strip():
            parts = line.split()
            if len(parts) >= 2:
                res = parts[1][:4].upper()
                if res:
                    return res
    raise ValueError(
        f"FATAL: Could not parse canonical residue name from @<TRIPOS>SUBSTRUCTURE in {mol2_file}. "
        "File is corrupted or lacks standard Tripos substructure records."
    )


def build_tleap_atom_types_block(mol2_files: list[str]) -> str:
    """
    Parses @<TRIPOS>ATOM and @<TRIPOS>BOND from custom OpenFF/GAFF2 mol2 files
    and constructs an explicit addAtomTypes declaration for tLEaP to prevent Sp0 hybridization warnings.
    """
    type_map: dict[str, tuple[str, str]] = {}
    for m_str in mol2_files:
        p = Path(m_str)
        if not p.exists():
            continue
        lines = p.read_text(encoding="utf-8", errors="ignore").splitlines()
        section = None
        atoms_info: dict[int, tuple[str, str]] = {}
        has_multiple_or_ar_bond: set[int] = set()

        for line in lines:
            l_strip = line.strip()
            if l_strip.startswith("@<TRIPOS>"):
                section = l_strip
                continue
            if not l_strip:
                continue
            if section == "@<TRIPOS>ATOM":
                parts = l_strip.split()
                if len(parts) >= 6 and parts[0].isdigit():
                    a_idx = int(parts[0])
                    a_name = parts[1]
                    a_type = parts[5]
                    clean_alpha = "".join([c for c in a_name if c.isalpha()])
                    if clean_alpha[:2].capitalize() in {"Cl", "Br", "Fe", "Mg", "Zn", "Ca", "Na"}:
                        elem = clean_alpha[:2].capitalize()
                    else:
                        elem = clean_alpha[:1].upper() if clean_alpha else "C"
                    atoms_info[a_idx] = (a_type, elem)
            elif section == "@<TRIPOS>BOND":
                parts = l_strip.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    a1, a2, b_type = int(parts[1]), int(parts[2]), parts[3].lower()
                    if b_type in {"2", "3", "ar", "am"}:
                        has_multiple_or_ar_bond.add(a1)
                        has_multiple_or_ar_bond.add(a2)

        for a_idx, (a_type, elem) in atoms_info.items():
            hyb = "sp2" if a_idx in has_multiple_or_ar_bond else "sp3"
            type_map[a_type] = (elem, hyb)

    if not type_map:
        return ""

    entries = [f'  {{ "{atype}" "{elem}" "{hyb}" }}' for atype, (elem, hyb) in sorted(type_map.items())]
    return "addAtomTypes {\n" + "\n".join(entries) + "\n}\n"


def cmd_antechamber(args):
    from rdkit import Chem

    res_name = args.res_name.strip().upper()[:4]
    if len(res_name) < 3:
        res_name = f"{res_name:<3s}"

    suppl = Chem.SDMolSupplier(str(args.sdf), removeHs=False, sanitize=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(f"FATAL: Could not read molecule from {args.sdf}")

    Chem.SanitizeMol(mol)

    clean_v2000_sdf = f"{args.mol_id}_v2000.sdf"
    with Chem.SDWriter(clean_v2000_sdf) as writer:
        writer.SetForceV3000(False)
        writer.write(mol)

    out_mol2 = f"{args.mol_id}_gaff2.mol2"
    scfconv = getattr(args, "scfconv", None) or "1.d-10"
    grms_tol = getattr(args, "grms_tol", None) or "0.0005"
    run_cmd(
        [
            "antechamber",
            "-i",
            clean_v2000_sdf,
            "-fi",
            "sdf",
            "-o",
            out_mol2,
            "-fo",
            "mol2",
            "-at",
            "gaff2",
            "-c",
            "bcc",
            "-nc",
            str(args.charge),
            "-rn",
            res_name,
            "-ek",
            f"qm_theory='AM1', grms_tol={grms_tol}, scfconv={scfconv}, ndiis_attempts=700",
        ]
    )
    run_cmd(
        [
            "parmchk2",
            "-i",
            out_mol2,
            "-f",
            "mol2",
            "-o",
            f"{args.mol_id}.frcmod",
            "-s",
            "2",
        ]
    )


def cmd_tleap_dry(args):
    """Generates pristine, un-solvated AMBER topology of the complex strictly for MM-PBSA."""
    subprocess.run(
        [
            sys.executable,
            f"{Path(__file__).parent}/amber_cap_fixer.py",
            args.complex_pdb,
            "complex_clean.pdb",
            "--bonds_out",
            "apply_bonds.in",
            "--ss_cutoff",
            str(args.ss_cutoff),
            "--water_names",
            "WAT,HOH,SOL,TIP3,OPC",
        ],
        check=True,
    )

    atom_types_block = build_tleap_atom_types_block(args.mol2s)
    with open("tleap_dry.in", "w") as f:
        f.write(
            f"source {args.protein_ff}\nsource {args.gaff_ff}\nsource {args.lipid_ff}\nsource {args.water_leaprc}\n"
        )
        if atom_types_block:
            f.write(atom_types_block)
        for m in args.mol2s:
            res_code = parse_mol2_resname(Path(m))
            f.write(f"{res_code} = loadmol2 {m}\n")
        for m in args.frcmods:
            f.write(f"loadAmberParams {m}\n")
        f.write(
            'addPdbResMap { {"Zn2+" "ZN"} {"Ca2+" "CA"} {"Mg2+" "MG"} {"Na+" "NA"} {"K+" "K"} {"Cl-" "CL"} }\n'
        )
        f.write('addPdbAtomMap { {"HZ" "HZ3"} }\n')
        f.write("complex = loadPdb complex_clean.pdb\n")
        f.write("set default PBRadii mbondi3\nsource apply_bonds.in\n")
        f.write("charge complex\n")
        f.write(
            "saveAmberParm complex complex_dry.prmtop complex_dry.inpcrd\nsavepdb complex complex_dry_tleap.pdb\nquit\n"
        )

    run_cmd(["tleap", "-f", "tleap_dry.in"])


def cmd_tleap(args):
    # For membrane systems, packmol-memgen has already solvated the bilayer with explicit water;
    # never strip packmol water in amber_cap_fixer!
    effective_water_strip = "" if args.membrane else args.water_names
    subprocess.run(
        [
            sys.executable,
            f"{Path(__file__).parent}/amber_cap_fixer.py",
            args.complex_pdb,
            "complex_clean.pdb",
            "--bonds_out",
            "apply_bonds.in",
            "--ss_cutoff",
            str(args.ss_cutoff),
            "--water_names",
            effective_water_strip,
        ],
        check=True,
    )

    box = args.water_box
    leaprc_water = args.water_leaprc
    box_cmd = (
        f"solvateBox complex {box} {args.box_padding}"
        if args.box_type == "box"
        else f"solvateOct complex {box} {args.box_padding}"
    )

    atom_types_block = build_tleap_atom_types_block(args.mol2s)

    # 1. First tLEaP pass: Solvation and accurate water molecule enumeration
    with open("tleap_solv.in", "w") as f:
        f.write(
            f"source {args.protein_ff}\nsource {args.gaff_ff}\nsource {args.lipid_ff}\nsource {leaprc_water}\n"
        )
        if atom_types_block:
            f.write(atom_types_block)
        for m in args.mol2s:
            res_code = parse_mol2_resname(Path(m))
            f.write(f"{res_code} = loadmol2 {m}\n")
        for m in args.frcmods:
            f.write(f"loadAmberParams {m}\n")
        f.write(
            'addPdbResMap { {"Zn2+" "ZN"} {"Ca2+" "CA"} {"Mg2+" "MG"} {"Na+" "NA"} {"K+" "K"} {"Cl-" "CL"} }\n'
        )
        f.write('addPdbAtomMap { {"HZ" "HZ3"} }\n')
        f.write("complex = loadPdb complex_clean.pdb\n")
        f.write("set default PBRadii mbondi3\nsource apply_bonds.in\n")
        if not args.membrane:
            f.write(f"{box_cmd}\n")
        else:
            f.write("setBox complex vdw\n")
        f.write("savepdb complex complex_solv.pdb\nquit\n")

    run_cmd(["tleap", "-f", "tleap_solv.in"])

    # 2. Compute exact physiological salt ion pairs from real water count
    salt_pairs = 0
    salt_conc = float(args.salt_conc) if getattr(args, "salt_conc", None) else 0.0
    if not args.membrane and salt_conc > 0.0 and Path("complex_solv.pdb").exists():
        water_residues = {
            w.strip().upper() for w in args.water_names.split(",") if w.strip()
        }
        wat_count = 0
        current_res_key = None
        for line in Path("complex_solv.pdb").read_text(errors="ignore").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                r_name = line[17:20].strip().upper()
                r_seq = line[21:26]
                if r_name in water_residues and r_seq != current_res_key:
                    current_res_key = r_seq
                    wat_count += 1
        WATER_MOLARITY = 55.5084  # Standard water molarity (mol/L) at 298.15 K
        if wat_count == 0:
            raise ValueError(
                "FATAL: Solvation failed! Real water count in complex_solv.pdb is zero."
            )
        salt_pairs = int(round(wat_count * (salt_conc / WATER_MOLARITY)))
        print(
            f"Adding physiological salt: {salt_conc} M ({salt_pairs} K+/Cl- pairs for {wat_count} water molecules)."
        )

    # 3. Final production tLEaP topology generation with neutralization + physiological salt
    with open("tleap.in", "w") as f:
        f.write(
            f"source {args.protein_ff}\nsource {args.gaff_ff}\nsource {args.lipid_ff}\nsource {leaprc_water}\n"
        )
        if atom_types_block:
            f.write(atom_types_block)
        for m in args.mol2s:
            m_path = Path(m)
            res_code = parse_mol2_resname(m_path)
            f.write(f"{res_code} = loadmol2 {m}\n")
        for m in args.frcmods:
            f.write(f"loadAmberParams {m}\n")
        f.write(
            'addPdbResMap { {"Zn2+" "ZN"} {"Ca2+" "CA"} {"Mg2+" "MG"} {"Na+" "NA"} {"K+" "K"} {"Cl-" "CL"} }\n'
        )
        f.write('addPdbAtomMap { {"HZ" "HZ3"} }\n')
        f.write("complex = loadPdb complex_clean.pdb\n")
        f.write("set default PBRadii mbondi3\nsource apply_bonds.in\n")
        if not args.membrane:
            f.write(f"{box_cmd}\n")
            f.write("addIonsRand complex Cl- 0\naddIonsRand complex K+ 0\n")
            if salt_pairs > 0:
                f.write(f"addIonsRand complex K+ {salt_pairs} Cl- {salt_pairs}\n")
        else:
            f.write("setBox complex vdw\n")
            f.write("addIonsRand complex Cl- 0\naddIonsRand complex K+ 0\n")
        f.write("charge complex\n")
        f.write(
            "saveAmberParm complex complex.prmtop complex.inpcrd\nsavepdb complex complex_tleap.pdb\nquit\n"
        )

    run_cmd(["tleap", "-f", "tleap.in"])


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    p_ante = subparsers.add_parser("antechamber")
    p_ante.add_argument("--sdf", required=True)
    p_ante.add_argument("--charge", required=True, type=int)
    p_ante.add_argument("--mol_id", required=True)
    p_ante.add_argument("--res_name", required=True)
    p_ante.add_argument("--scfconv", default="1.d-10")
    p_ante.add_argument("--grms_tol", default="0.0005")

    p_tleap = subparsers.add_parser("tleap")
    p_tleap.add_argument("--complex_pdb", required=True)
    p_tleap.add_argument("--water_leaprc", required=True)
    p_tleap.add_argument("--water_box", required=True)
    p_tleap.add_argument("--box_type", required=True)
    p_tleap.add_argument("--box_padding", required=True)
    p_tleap.add_argument("--ss_cutoff", required=True)
    p_tleap.add_argument("--protein_ff", required=True)
    p_tleap.add_argument("--gaff_ff", required=True)
    p_tleap.add_argument("--lipid_ff", required=True)
    p_tleap.add_argument("--water_names", required=True)
    p_tleap.add_argument("--salt_conc", required=True)
    p_tleap.add_argument("--mol2s", nargs="*", default=[])
    p_tleap.add_argument("--frcmods", nargs="*", default=[])
    p_tleap.add_argument("--membrane", action="store_true")

    p_dry = subparsers.add_parser("tleap_dry")
    p_dry.add_argument("--complex_pdb", required=True)
    p_dry.add_argument("--water_leaprc", required=True)
    p_dry.add_argument("--ss_cutoff", required=True)
    p_dry.add_argument("--protein_ff", required=True)
    p_dry.add_argument("--gaff_ff", required=True)
    p_dry.add_argument("--lipid_ff", required=True)
    p_dry.add_argument("--mol2s", nargs="*", default=[])
    p_dry.add_argument("--frcmods", nargs="*", default=[])

    args = parser.parse_args()
    if args.subcommand == "antechamber":
        cmd_antechamber(args)
    elif args.subcommand == "tleap":
        cmd_tleap(args)
    elif args.subcommand == "tleap_dry":
        cmd_tleap_dry(args)


if __name__ == "__main__":
    main()
