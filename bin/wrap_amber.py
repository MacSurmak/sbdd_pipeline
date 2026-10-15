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


def cmd_antechamber(args):
    res_name = args.res_name.strip().upper()[:4]
    if len(res_name) < 3:
        res_name = f"{res_name:<3s}"
    
    out_mol2 = f"{args.mol_id}_gaff2.mol2"
    run_cmd(
        [
            "antechamber",
            "-i", args.sdf,
            "-fi", "sdf",
            "-o", out_mol2,
            "-fo", "mol2",
            "-at", "gaff2",
            "-c", "bcc",
            "-nc", str(args.charge),
            "-rn", res_name,
            "-ek", "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, MAXCYC=0",
        ]
    )
    run_cmd(
        [
            "parmchk2",
            "-i", f"{args.mol_id}_gaff2.mol2",
            "-f", "mol2",
            "-o", f"{args.mol_id}.frcmod",
            "-s", "2",
        ]
    )


def cmd_tleap(args):
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
            args.water_names,
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

    # 1. First tLEaP pass: Solvation and accurate water molecule enumeration
    with open("tleap_solv.in", "w") as f:
        f.write(
            f"source leaprc.protein.ff19SB\nsource leaprc.gaff2\nsource leaprc.lipid21\nsource {leaprc_water}\n"
        )
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
            f"source leaprc.protein.ff19SB\nsource leaprc.gaff2\nsource leaprc.lipid21\nsource {leaprc_water}\n"
        )
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

    p_tleap = subparsers.add_parser("tleap")
    p_tleap.add_argument("--complex_pdb", required=True)
    p_tleap.add_argument("--water_leaprc", required=True)
    p_tleap.add_argument("--water_box", required=True)
    p_tleap.add_argument("--box_type", required=True)
    p_tleap.add_argument("--box_padding", required=True)
    p_tleap.add_argument("--ss_cutoff", required=True)
    p_tleap.add_argument("--water_names", required=True)
    p_tleap.add_argument("--salt_conc", default="0.15")
    p_tleap.add_argument("--mol2s", nargs="*", default=[])
    p_tleap.add_argument("--frcmods", nargs="*", default=[])
    p_tleap.add_argument("--membrane", action="store_true")

    args = parser.parse_args()
    if args.subcommand == "antechamber":
        cmd_antechamber(args)
    elif args.subcommand == "tleap":
        cmd_tleap(args)


if __name__ == "__main__":
    main()
