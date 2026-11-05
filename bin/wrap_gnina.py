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

    import os
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")

    templates_by_idx: dict[int, tuple[str, Chem.Mol]] = {}
    effective_ligand_input = args.ligand
    tagged_input_sdf = f"{args.out_sdf}.in_tagged.sdf"

    if args.mode == "dock":
        with Chem.SDWriter(tagged_input_sdf) as tag_w:
            tag_w.SetForceV3000(False)
            for idx, tmpl_mol in enumerate(
                Chem.SDMolSupplier(args.ligand, removeHs=False, sanitize=False)
            ):
                if tmpl_mol is None:
                    continue
                t_name = (
                    tmpl_mol.GetProp("_Name").strip()
                    if tmpl_mol.HasProp("_Name")
                    else ""
                )
                if not t_name and tmpl_mol.HasProp("s_m_title"):
                    t_name = tmpl_mol.GetProp("s_m_title").strip()
                if not t_name:
                    raise ValueError(
                        f"FATAL: Unnamed molecule at index {idx} in {args.ligand}!"
                    )
                tmpl_mol.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(tmpl_mol)
                templates_by_idx[idx] = (t_name, Chem.Mol(tmpl_mol))

                tagged_mol = Chem.Mol(tmpl_mol)
                tagged_name = f"{t_name}__tmpl_{idx}"
                tagged_mol.SetProp("_Name", tagged_name)
                if tagged_mol.HasProp("s_m_title"):
                    tagged_mol.SetProp("s_m_title", tagged_name)
                tag_w.write(tagged_mol)
        effective_ligand_input = tagged_input_sdf

    cmd = [
        "gnina",
        "-r",
        args.receptor,
        "-l",
        effective_ligand_input,
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
        Path(tagged_input_sdf).unlink(missing_ok=True)
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
        framework_dir = Path(__file__).resolve().parent.parent
        if str(framework_dir) not in sys.path:
            sys.path.insert(0, str(framework_dir))
        from cli.core.chem import transfer_coordinates_to_template

        out_suppl = Chem.SDMolSupplier(args.out_sdf, removeHs=False, sanitize=False)
        tmp_sdf = f"{args.out_sdf}.tmp"

        if args.mode == "rescore":
            in_suppl = Chem.SDMolSupplier(args.ligand, removeHs=False, sanitize=False)
            with Chem.SDWriter(tmp_sdf) as w:
                w.SetForceV3000(False)
                for orig_m, sc_m in zip(in_suppl, out_suppl):
                    if orig_m is None or sc_m is None:
                        continue
                    for prop in sc_m.GetPropNames():
                        orig_m.SetProp(prop, sc_m.GetProp(prop))
                    w.write(orig_m)
        else:
            import json
            import socket
            import time

            raw_poses = [m for m in out_suppl if m is not None]
            total_poses = len(raw_poses)
            udp_sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            master_ip = os.environ.get("MASTER_IP", "192.168.1.101")
            udp_port = int(os.environ.get("UDP_PORT", "8089"))
            node_name = (
                os.environ.get("SLURM_NODENAME")
                or os.environ.get("SLURMD_NODENAME")
                or socket.gethostname().split(".")[0]
            )
            task_id = (
                os.environ.get("NXF_TASK_NAME")
                or (
                    f"{os.environ.get('SLURM_JOB_NAME')} (#{os.environ.get('SLURM_JOB_ID', 0)})"
                    if os.environ.get("SLURM_JOB_NAME")
                    else None
                )
                or f"gnina_post_{os.getpid()}"
            )
            t_post_start = time.monotonic()

            def _emit_post_udp(p_idx: int, ev: str = "update"):
                if total_poses == 0:
                    return
                el = max(0.001, time.monotonic() - t_post_start)
                rate = p_idx / el if p_idx > 0 else 0.0
                rem = max(0, total_poses - p_idx)
                eta_s = f"{int(rem / rate)}s" if rate > 0 and rem > 0 else "done"
                payload = {
                    "tool": "gnina",
                    "task": task_id,
                    "session_id": f"gnina_{task_id}",
                    "stage": "Restore H",
                    "node": node_name,
                    "event": ev,
                    "step": p_idx,
                    "total": total_poses,
                    "dt_ps": None,
                    "unit": "pose",
                    "hardware": {"cpu": {"allocated_cores": args.cpus}},
                    "performance": {"elapsed_s": round(el, 1)},
                    "display_step": p_idx,
                    "display_total": total_poses,
                    "step_info": f"[cyan]{p_idx}/{total_poses} pose[/cyan] [dim](Restore H)[/dim]",
                    "speed": f"{rate:>5.1f} pose/s" if rate > 0 else "-- pose/s",
                    "eta": eta_s,
                }
                try:
                    udp_sock.sendto(
                        json.dumps(payload).encode("utf-8"), (master_ip, udp_port)
                    )
                except Exception:
                    pass

            with Chem.SDWriter(tmp_sdf) as w:
                w.SetForceV3000(False)
                for p_idx, m in enumerate(raw_poses, 1):
                    raw_tag = m.GetProp("_Name").strip() if m.HasProp("_Name") else ""
                    if not raw_tag and m.HasProp("s_m_title"):
                        raw_tag = m.GetProp("s_m_title").strip()
                    if "__tmpl_" not in raw_tag:
                        raise ValueError(
                            f"FATAL: Missing template index tag in GNINA output pose '{raw_tag}'!"
                        )
                    orig_title, idx_str = raw_tag.rsplit("__tmpl_", 1)
                    tmpl_entry = templates_by_idx.get(int(idx_str))
                    if tmpl_entry is None:
                        raise KeyError(
                            f"FATAL: Template index {idx_str} for '{orig_title}' not found in {args.ligand}!"
                        )
                    _, exact_tmpl = tmpl_entry
                    restored_mol = transfer_coordinates_to_template(
                        template_mol=exact_tmpl,
                        source_mol=m,
                        relax_missing_hydrogens=True,
                    )
                    for prop in m.GetPropNames():
                        restored_mol.SetProp(prop, m.GetProp(prop))
                    restored_mol.SetProp("_Name", orig_title)
                    if restored_mol.HasProp("s_m_title"):
                        restored_mol.SetProp("s_m_title", orig_title)
                    w.write(restored_mol)
                    if p_idx % 5 == 0 or p_idx == total_poses:
                        _emit_post_udp(p_idx)

            _emit_post_udp(total_poses, ev="exit")

        os.replace(tmp_sdf, args.out_sdf)
        Path(tagged_input_sdf).unlink(missing_ok=True)
    except Exception as e:
        print(
            f"FATAL: Failed to restore and normalize GNINA SDF {args.out_sdf} via RDKit: {e}",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
