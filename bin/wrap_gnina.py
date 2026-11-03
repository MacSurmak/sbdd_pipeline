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
        if args.cnn_model:
            cmd.extend(["--cnn", args.cnn_model])
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

    # Restore full-atom chemical topology and explicit hydrogens directly in worker $SLURM_TMPDIR
    try:
        import os
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        framework_dir = Path(__file__).resolve().parent.parent
        if str(framework_dir) not in sys.path:
            sys.path.insert(0, str(framework_dir))
        from cli.core.chem import transfer_coordinates_to_template

        out_suppl = Chem.SDMolSupplier(args.out_sdf, removeHs=False, sanitize=False)
        tmp_sdf = f"{args.out_sdf}.tmp"

        if args.mode == "rescore":
            in_suppl = Chem.SDMolSupplier(args.ligand, removeHs=False, sanitize=False)
            in_mols = [m for m in in_suppl if m is not None]
            scored_mols = [m for m in out_suppl if m is not None]
            if len(in_mols) != len(scored_mols):
                raise ValueError(
                    f"FATAL: GNINA rescore count mismatch: input had {len(in_mols)} poses, output has {len(scored_mols)}"
                )
            with Chem.SDWriter(tmp_sdf) as w:
                w.SetForceV3000(False)
                for orig_m, sc_m in zip(in_mols, scored_mols):
                    for prop in sc_m.GetPropNames():
                        orig_m.SetProp(prop, sc_m.GetProp(prop))
                    w.write(orig_m)
        else:
            templates_by_title: dict[str, list[Chem.Mol]] = {}
            for tmpl_mol in Chem.SDMolSupplier(args.ligand, removeHs=False, sanitize=False):
                if tmpl_mol is not None:
                    t_name = tmpl_mol.GetProp("_Name").strip() if tmpl_mol.HasProp("_Name") else ""
                    if not t_name and tmpl_mol.HasProp("s_m_title"):
                        t_name = tmpl_mol.GetProp("s_m_title").strip()
                    if t_name:
                        tmpl_mol.UpdatePropertyCache(strict=False)
                        Chem.SanitizeMol(tmpl_mol)
                        templates_by_title.setdefault(t_name, []).append(tmpl_mol)

            with Chem.SDWriter(tmp_sdf) as w:
                w.SetForceV3000(False)
                for m in out_suppl:
                    if m is None:
                        continue
                    m_title = m.GetProp("_Name").strip() if m.HasProp("_Name") else ""
                    if not m_title and m.HasProp("s_m_title"):
                        m_title = m.GetProp("s_m_title").strip()
                    candidates = templates_by_title.get(m_title)
                    if not candidates:
                        raise KeyError(
                            f"FATAL: Docked pose title '{m_title}' not found in input chunk {args.ligand}!"
                        )
                    restored_mol = None
                    last_err = None
                    for cand_tmpl in candidates:
                        try:
                            restored_mol = transfer_coordinates_to_template(
                                template_mol=cand_tmpl,
                                source_mol=m,
                                relax_missing_hydrogens=True,
                            )
                            break
                        except ValueError as err:
                            last_err = err
                            continue
                    if restored_mol is None:
                        raise ValueError(
                            f"FATAL: Failed to restore all-atom topology for GNINA pose '{m_title}': {last_err}"
                        )
                    w.write(restored_mol)

        os.replace(tmp_sdf, args.out_sdf)
    except Exception as e:
        print(
            f"FATAL: Failed to restore and normalize GNINA SDF {args.out_sdf} via RDKit: {e}",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
