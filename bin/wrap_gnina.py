#!/usr/bin/env python3
import argparse
import itertools
import os
import subprocess
import sys
from pathlib import Path
from typing import Final
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

FRAMEWORK_DIR: Final[Path] = Path(__file__).resolve().parent.parent
if str(FRAMEWORK_DIR) not in sys.path:
    sys.path.insert(0, str(FRAMEWORK_DIR))

from cli.core.chem import transfer_coordinates_to_template


def main() -> None:
    parser = argparse.ArgumentParser(description="Universal GNINA 3D-CNN Pose Rescorer")
    parser.add_argument("mode", choices=["rescore"])
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--ref_ligand", required=True)
    parser.add_argument("--cpus", type=int, required=True)
    parser.add_argument("--rescore_mode", required=True, choices=["score_only", "minimize"])
    parser.add_argument("--autobox_add", required=True, type=float)
    parser.add_argument("--cnn_model", required=True)
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
        "--autobox_add",
        str(args.autobox_add),
        "--cpu",
        str(args.cpus),
    ]

    if args.rescore_mode == "minimize":
        cmd.extend(["--minimize", "--cnn_scoring", "rescore"])
    else:
        cmd.extend(["--score_only", "--cnn_scoring", "rescore"])

    if args.cnn_model and args.cnn_model != "default":
        cmd.extend(["--cnn", args.cnn_model])

    cmd.extend(["-o", args.out_sdf])

    print(f"Exec: {' '.join(cmd)}", flush=True)
    with open(args.out_log, "w", encoding="utf-8") as f:
        res = subprocess.run(
            cmd,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=False,
            env=os.environ.copy(),
            timeout=14400,
        )

    if res.returncode != 0:
        log_content = ""
        if Path(args.out_log).exists():
            log_content = Path(args.out_log).read_text(errors="ignore")
        print(
            f"FATAL: GNINA rescoring failed with exit code {res.returncode}!\n"
            f"Command: {' '.join(cmd)}\n"
            f"==================== GNINA LOG ({args.out_log}) ====================\n"
            f"{log_content}\n"
            f"====================================================================",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(res.returncode)

    in_suppl = Chem.SDMolSupplier(args.ligand, removeHs=False, sanitize=False)
    out_suppl = Chem.SDMolSupplier(args.out_sdf, removeHs=False, sanitize=False)

    tmp_sdf = f"{args.out_sdf}.tmp"
    pose_count = 0
    with Chem.SDWriter(tmp_sdf) as w:
        w.SetForceV3000(False)
        for idx, (orig_m, sc_m) in enumerate(
            itertools.zip_longest(in_suppl, out_suppl, fillvalue=None), start=1
        ):
            if orig_m is None or sc_m is None:
                raise ValueError(
                    f"FATAL: Pose count mismatch or corrupted molecule record at pose index {idx} "
                    f"during GNINA rescoring of {args.ligand}!"
                )
            if args.rescore_mode == "minimize":
                orig_m.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(orig_m)
                final_m = transfer_coordinates_to_template(
                    template_mol=orig_m,
                    source_mol=sc_m,
                    relax_missing_hydrogens=True,
                )
                for p_name in orig_m.GetPropNames():
                    final_m.SetProp(p_name, orig_m.GetProp(p_name))
            else:
                final_m = orig_m

            orig_unidock_score = (
                final_m.GetProp("unidock_vina_score")
                if final_m.HasProp("unidock_vina_score")
                else None
            )
            for prop in sc_m.GetPropNames():
                final_m.SetProp(prop, sc_m.GetProp(prop))
            if orig_unidock_score is not None:
                final_m.SetProp("unidock_vina_score", orig_unidock_score)
            w.write(final_m)
            pose_count += 1

    if pose_count == 0:
        raise ValueError(
            f"FATAL: Zero poses processed during GNINA rescoring of {args.ligand}!"
        )

    os.replace(tmp_sdf, args.out_sdf)


if __name__ == "__main__":
    main()
