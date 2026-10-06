#!/usr/bin/env python3
import argparse
import sys
import subprocess
import shutil
from pathlib import Path


def run_cmd(cmd, stdin=None, stdout=None):
    print(f"Exec: {' '.join(cmd)}")
    res = subprocess.run(
        cmd, stdin=stdin, stdout=stdout, stderr=subprocess.STDOUT, check=False
    )
    if res.returncode != 0:
        print(f"FATAL: Command failed: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(res.returncode)


def cmd_packmol(args):
    cmd = [
        "packmol-memgen",
        "--pdb",
        args.pdb,
        "--lipids",
        f"{args.upper_lipids}//{args.lower_lipids}",
        "--ratio",
        f"{args.upper_ratio}//{args.lower_ratio}",
        "--preoriented",
        "--dist",
        args.lipid_pad,
        "--dist_wat",
        args.wat_pad,
        "--salt",
        "--salt_c",
        "K+",
        "--saltcon",
        args.salt_conc,
        "--notprotonate",
        "-o",
        f"{args.out_id}_mem_complex",
    ]
    run_cmd(cmd)

    pdb1 = Path(f"{args.out_id}_mem_complex.pdb")
    pdb2 = Path(f"{args.out_id}_mem_complex_mempack.pdb")
    if pdb2.exists():
        pdb2.rename(pdb1)
    if not pdb1.exists():
        print("FATAL: Packmol failed to generate PDB", file=sys.stderr)
        sys.exit(1)


def cmd_immers(args):
    res_lib = Path("/usr/local/share/immers/res.lib")
    if res_lib.exists():
        Path("res.lib").symlink_to(res_lib)

    with open("immers.inp", "w") as f:
        f.write(
            f"2\nno\n{args.pdb}\n1\n{args.membrane_type}\nplanar\n{args.membrane_topo}\n{args.membrane_chains}\n"
        )

    immers_bin = (
        shutil.which("immers")
        or "/opt/containers/bin/immers"
        or "/usr/local/bin/immers"
    )
    with open("immers.inp", "r") as fin, open("immers.log", "w") as fout:
        run_cmd([immers_bin], stdin=fin, stdout=fout)

    pdb_base = Path(args.pdb).stem
    out1, out2, final = (
        Path(f"{pdb_base}out.pdb"),
        Path("protout.pdb"),
        Path(f"{args.out_id}_oriented.pdb"),
    )

    if out1.exists():
        out1.rename(final)
    elif out2.exists():
        out2.rename(final)
    else:
        print("FATAL: Immers PDB not found", file=sys.stderr)
        sys.exit(1)

    subprocess.run(
        [
            sys.executable,
            f"{Path(__file__).parent}/parse_immers.py",
            "--datapar",
            "datapar1",
            "--out_json",
            f"{args.out_id}_membrane.json",
            "--out_txt",
            f"{args.out_id}_membrane.txt",
        ],
        check=True,
    )


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    p_pack = subparsers.add_parser("packmol")
    p_pack.add_argument("--pdb", required=True)
    p_pack.add_argument("--upper_lipids", required=True)
    p_pack.add_argument("--upper_ratio", required=True)
    p_pack.add_argument("--lower_lipids", required=True)
    p_pack.add_argument("--lower_ratio", required=True)
    p_pack.add_argument("--lipid_pad", required=True)
    p_pack.add_argument("--wat_pad", required=True)
    p_pack.add_argument("--salt_conc", required=True)
    p_pack.add_argument("--out_id", required=True)

    p_imm = subparsers.add_parser("immers")
    p_imm.add_argument("--pdb", required=True)
    p_imm.add_argument("--membrane_type", required=True)
    p_imm.add_argument("--membrane_topo", required=True)
    p_imm.add_argument("--membrane_chains", required=True)
    p_imm.add_argument("--out_id", required=True)

    args = parser.parse_args()
    if args.subcommand == "packmol":
        cmd_packmol(args)
    elif args.subcommand == "immers":
        cmd_immers(args)


if __name__ == "__main__":
    main()
