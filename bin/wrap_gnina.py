#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["dock", "rescore"])
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--ref_ligand", required=True)
    parser.add_argument("--cpus", type=int, required=True)
    parser.add_argument("--autobox_add")
    parser.add_argument("--exhaustiveness")
    parser.add_argument("--cnn_scoring")
    parser.add_argument("--num_modes")
    parser.add_argument("--cnn_model")
    parser.add_argument("--out_sdf", required=True)
    parser.add_argument("--out_log", required=True)
    args = parser.parse_args()

    ref = Path(args.ref_ligand)
    if not ref.exists() or ref.stat().st_size == 0:
        print(
            f"FATAL: Reference ligand missing or empty: {args.ref_ligand}",
            file=sys.stderr,
        )
        sys.exit(1)

    cmd = [
        "gnina",
        "-r",
        args.receptor,
        "-l",
        args.ligand,
        "--autobox_ligand",
        args.ref_ligand,
        "--cpu",
        str(args.cpus),
    ]

    if args.mode == "dock":
        cmd.extend(
            [
                "--autobox_add",
                args.autobox_add,
                "--exhaustiveness",
                args.exhaustiveness,
                "--cnn_scoring",
                args.cnn_scoring,
                "--num_modes",
                args.num_modes,
            ]
        )
    elif args.mode == "rescore":
        cmd.extend(["--score_only"])
        if args.cnn_model:
            cmd.extend(["--cnn", args.cnn_model])

    cmd.extend(["-o", args.out_sdf])

    print(f"Exec: {' '.join(cmd)}", flush=True)
    with open(args.out_log, "w", encoding="utf-8") as f:
        res = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=False)

    if res.returncode != 0:
        log_content = ""
        if Path(args.out_log).exists():
            log_content = Path(args.out_log).read_text(errors="ignore")
        print(
            f"FATAL: GNINA failed with exit code {res.returncode}!\n"
            f"Command: {' '.join(cmd)}\n"
            f"==================== GNINA LOG ({args.out_log}) ====================\n"
            f"{log_content}\n"
            f"====================================================================",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(res.returncode)


if __name__ == "__main__":
    main()
