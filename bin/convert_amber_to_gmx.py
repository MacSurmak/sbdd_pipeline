#!/usr/bin/env python3
"""
bin/convert_amber_to_gmx.py
Converts Amber PRMTOP and INPCRD files to GROMACS topology and GRO file using ParmEd.
"""
import argparse
import sys
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Convert Amber Topology to GROMACS via ParmEd"
    )
    parser.add_argument("--prmtop", required=True, help="Input Amber prmtop file")
    parser.add_argument("--inpcrd", default=None, help="Input Amber inpcrd file")
    parser.add_argument("--out_top", default="topol.top", help="Output GROMACS top file")
    parser.add_argument("--out_gro", default="step3_input.gro", help="Output GROMACS gro file")
    parser.add_argument("--hmr", action="store_true", help="Enable Hydrogen Mass Repartitioning")
    parser.add_argument("--hmr_mass", type=float, default=None, help="Target hydrogen mass in Da for HMR")
    args = parser.parse_args()

    try:
        import parmed as pmd
    except ImportError:
        print(
            "Error: ParmEd is not installed in the current environment.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.inpcrd:
        print(f"Loading {args.prmtop} and {args.inpcrd}...")
        struct = pmd.load_file(args.prmtop, args.inpcrd)
    else:
        print(f"Loading dry topology {args.prmtop}...")
        struct = pmd.load_file(args.prmtop)

    # Sanitize atomic_number on both Atom and AtomType objects strictly from initial physical masses
    # in prmtop so both ParmEd HMR and GROMACS [ atomtypes ] column 2 receive exact element numbers
    for a in struct.atoms:
        m = float(a.mass)
        if m < 0.5:
            at_num = 0
        elif 0.8 <= m <= 2.5:
            at_num = 1
        elif m < 13.0:
            at_num = 6
        elif m < 15.0:
            at_num = 7
        elif m < 17.5:
            at_num = 8
        elif m < 20.0:
            at_num = 9
        elif m < 23.5:
            at_num = 11
        elif m < 25.5:
            at_num = 12
        elif m < 31.5:
            at_num = 15
        elif m < 33.5:
            at_num = 16
        elif m < 37.0:
            at_num = 17
        elif m < 41.0:
            at_num = 19
        elif m < 85.0:
            at_num = 35
        else:
            at_num = 53
        a.atomic_number = at_num
        if a.atom_type is not None:
            a.atom_type.atomic_number = at_num

    if args.inpcrd:
        if args.hmr:
            if args.hmr_mass is None or args.hmr_mass <= 1.008:
                print(
                    f"FATAL: --hmr requires a valid --hmr_mass > 1.008 Da, got: {args.hmr_mass}",
                    file=sys.stderr,
                )
                sys.exit(1)
            from parmed.tools import HMassRepartition

            print(
                f"Applying Hydrogen Mass Repartitioning (HMR = {args.hmr_mass:.3f} Da, dowater=False)..."
            )
            HMassRepartition(struct, args.hmr_mass).execute()

            water_resnames = {"WAT", "HOH", "SOL", "OPC", "OPC3", "TIP3", "TP3"}
            for a in struct.atoms:
                res_name = a.residue.name.strip().upper() if a.residue else ""
                if res_name in water_resnames or a.atomic_number == 0:
                    continue
                if a.atomic_number == 1 and abs(float(a.mass) - args.hmr_mass) > 0.05:
                    print(
                        f"FATAL: HMR verification failed on hydrogen atom {a.name} (res {res_name}): "
                        f"expected mass {args.hmr_mass:.3f} Da, got {a.mass:.3f} Da!",
                        file=sys.stderr,
                    )
                    sys.exit(1)
                if a.atomic_number > 1 and float(a.mass) < 4.5:
                    print(
                        f"FATAL: HMR verification failed on heavy atom {a.name} (res {res_name}): "
                        f"unphysical mass {a.mass:.3f} Da < 4.5 Da!",
                        file=sys.stderr,
                    )
                    sys.exit(1)

        # Strict Fail-Fast charge neutrality verification and 6-decimal rounding dust elimination for PME
        total_q = sum(float(a.charge) for a in struct.atoms)
        if abs(total_q) > 0.05:
            print(
                f"FATAL: Solvated system has non-zero net charge ({total_q:+.6f} e)! "
                "Ewald/PME electrostatics strictly require a neutral box. Execution halted by Fail-Fast law.",
                file=sys.stderr,
            )
            sys.exit(1)
        for a in struct.atoms:
            a.charge = round(float(a.charge), 6)
        rem_q = round(sum(a.charge for a in struct.atoms), 6)
        if abs(rem_q) > 0.0 and len(struct.atoms) > 0:
            struct.atoms[0].charge = round(struct.atoms[0].charge - rem_q, 6)

    print(f"Saving GROMACS {args.out_top}...")
    struct.save(args.out_top, format="gromacs", overwrite=True)

    # Fail-Fast Validation: Ensure terminal caps are covalently integrated and did not fragment
    topol_text = Path(args.out_top).read_text(encoding="utf-8", errors="ignore")
    isolated_caps = []
    in_moleculetype = False

    for raw_line in topol_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("["):
            in_moleculetype = "moleculetype" in line.lower()
            continue
        if in_moleculetype:
            if line.startswith(";") or line.startswith("#"):
                continue
            mname = line.split()[0].strip()
            if mname in ["NME", "NMA", "ACE"]:
                isolated_caps.append(mname)
            in_moleculetype = False

    if isolated_caps:
        print(
            f"FATAL: TopologyIntegrityError: Terminal cap(s) {isolated_caps} were exported as detached, "
            f"independent moleculetype(s) in topol.top! Peptide bond to protein chain was not formed in tLEaP. "
            f"Execution halted by Fail-Fast law.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.inpcrd:
        if struct.box is None or any(float(x) <= 0.0 for x in struct.box[:3]):
            print(
                f"FATAL: Periodic box vectors are missing or non-positive in {args.prmtop} / {args.inpcrd}: {struct.box}! "
                "Execution halted by Fail-Fast law.",
                file=sys.stderr,
            )
            sys.exit(1)
        print(
            f"Saving GROMACS {args.out_gro} (box: {struct.box[0]:.2f} x {struct.box[1]:.2f} x {struct.box[2]:.2f} A)..."
        )
        struct.save(args.out_gro, overwrite=True)
    print("Conversion successful.")


if __name__ == "__main__":
    main()
