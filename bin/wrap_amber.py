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

def cmd_antechamber(args):
    with open(args.charge_file, "r") as f:
        charge = f.read().strip()
    run_cmd(["antechamber", "-i", args.mol2, "-fi", "mol2", "-o", f"{args.mol_id}_gaff2.mol2", "-fo", "mol2", "-at", "gaff2", "-c", "bcc", "-nc", charge, "-rn", args.mol_id])
    run_cmd(["antechamber", "-i", f"{args.mol_id}_gaff2.mol2", "-fi", "mol2", "-o", f"{args.mol_id}_gaff2.pdb", "-fo", "pdb"])
    run_cmd(["parmchk2", "-i", f"{args.mol_id}_gaff2.mol2", "-f", "mol2", "-o", f"{args.mol_id}.frcmod", "-s", "2"])

def cmd_tleap(args):
    subprocess.run([sys.executable, f"{Path(__file__).parent}/amber_cap_fixer.py", args.complex_pdb, "complex_clean.pdb", "--bonds_out", "apply_bonds.in"], check=True)
    with open("tleap.in", "w") as f:
        f.write("source leaprc.protein.ff19SB\nsource leaprc.gaff2\nsource leaprc.lipid21\n")
        if args.water_model.lower() == "tip3p":
            f.write("source leaprc.water.tip3p\nloadAmberParams frcmod.ionslm_126_tip3p\n")
            box = "TIP3PBOX"
        else:
            f.write("source leaprc.water.opc\nloadAmberParams frcmod.ionslm_126_opc\n")
            box = "OPCBOX"
            
        for m in args.mol2s:
            name = Path(m).stem.replace("_gaff2", "")
            f.write(f"{name} = loadmol2 {m}\n")
        for m in args.frcmods:
            f.write(f"loadAmberParams {m}\n")
            
        f.write('addPdbResMap { {"Zn2+" "ZN"} {"Ca2+" "CA"} {"Mg2+" "MG"} {"Na+" "NA"} {"K+" "K"} {"Cl-" "CL"} }\n')
        f.write('addPdbAtomMap { {"HZ" "HZ3"} }\n')
        f.write("complex = loadPdb complex_clean.pdb\n")
        f.write("set default PBRadii mbondi3\nsource apply_bonds.in\n")
        if not args.membrane:
            f.write(f"solvateOct complex {box} {args.box_padding}\n")
        f.write("addIonsRand complex Cl- 0\naddIonsRand complex K+ 0\ncharge complex\n")
        f.write("saveAmberParm complex complex.prmtop complex.inpcrd\nsavepdb complex complex_tleap.pdb\nquit\n")
        
    run_cmd(["tleap", "-f", "tleap.in"])

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    p_ante = subparsers.add_parser("antechamber")
    p_ante.add_argument("--mol2", required=True)
    p_ante.add_argument("--charge_file", required=True)
    p_ante.add_argument("--mol_id", required=True)

    p_tleap = subparsers.add_parser("tleap")
    p_tleap.add_argument("--complex_pdb", required=True)
    p_tleap.add_argument("--water_model", required=True)
    p_tleap.add_argument("--box_padding", required=True)
    p_tleap.add_argument("--mol2s", nargs="*", default=[])
    p_tleap.add_argument("--frcmods", nargs="*", default=[])
    p_tleap.add_argument("--membrane", action="store_true")

    args = parser.parse_args()
    if args.subcommand == "antechamber": cmd_antechamber(args)
    elif args.subcommand == "tleap": cmd_tleap(args)

if __name__ == "__main__":
    main()
