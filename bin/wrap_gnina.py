#!/usr/bin/env python3
"""
bin/wrap_gnina.py
Executes GNINA docking or rescoring with explicit CPU thread control (--cpu)
and intra-batch parallel worker pool sharing the allocated GPU.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import sys
from pathlib import Path


def split_sdf_records(sdf_path: Path) -> list[str]:
    records = []
    current_lines = []
    with open(sdf_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.strip() == "$$$$":
                if any(l.strip() for l in current_lines):
                    records.append("".join(current_lines) + "$$$$\n")
                current_lines = []
            else:
                current_lines.append(line)
    if any(l.strip() for l in current_lines):
        records.append("".join(current_lines) + "\n$$$$\n")
    return records


def run_single_gnina_worker(
    worker_idx: int,
    in_sdf: Path,
    out_sdf: Path,
    out_log: Path,
    args,
    threads_per_worker: int,
) -> tuple[int, Path, Path]:
    cmd = [
        "gnina",
        "-r",
        args.receptor,
        "-l",
        str(in_sdf),
        "--autobox_ligand",
        args.ref_ligand,
        "--cpu",
        str(threads_per_worker),
    ]

    if args.mode == "dock":
        cmd.extend(
            [
                "--autobox_add",
                str(args.autobox_add),
                "--exhaustiveness",
                str(args.exhaustiveness),
                "--cnn_scoring",
                str(args.cnn_scoring),
                "--num_modes",
                str(args.num_modes),
            ]
        )
    elif args.mode == "rescore":
        cmd.extend(["--score_only"])
        if args.cnn_model:
            cmd.extend(["--cnn", str(args.cnn_model)])

    cmd.extend(["-o", str(out_sdf)])

    print(f"[Worker {worker_idx}] Exec: {' '.join(cmd)}", flush=True)
    with open(out_log, "w", encoding="utf-8") as f:
        res = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=False)

    if res.returncode != 0:
        log_content = (
            out_log.read_text(encoding="utf-8", errors="ignore")
            if out_log.exists()
            else ""
        )
        raise RuntimeError(
            f"FATAL: GNINA worker {worker_idx} failed with exit code {res.returncode}!\n"
            f"Command: {' '.join(cmd)}\n"
            f"==================== GNINA LOG ({out_log.name}) ====================\n"
            f"{log_content}\n"
            f"===================================================================="
        )
    if not out_sdf.exists() or out_sdf.stat().st_size == 0:
        raise RuntimeError(
            f"FATAL: GNINA worker {worker_idx} exited 0 but failed to produce {out_sdf}!"
        )
    return worker_idx, out_sdf, out_log


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["dock", "rescore"])
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--ref_ligand", required=True)
    parser.add_argument("--cpus", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
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

    lig_path = Path(args.ligand)
    if not lig_path.exists() or lig_path.stat().st_size == 0:
        print(f"FATAL: Input ligand SDF missing or empty: {lig_path}", file=sys.stderr)
        sys.exit(1)

    records = split_sdf_records(lig_path)
    if not records:
        print(f"FATAL: Zero molecule records parsed from {lig_path}", file=sys.stderr)
        sys.exit(1)

    n_workers = max(1, min(int(args.workers), len(records)))
    threads_per_worker = max(1, int(args.cpus) // n_workers)

    if n_workers == 1:
        try:
            run_single_gnina_worker(
                0,
                lig_path,
                Path(args.out_sdf),
                Path(args.out_log),
                args,
                threads_per_worker,
            )
        except Exception as e:
            print(str(e), file=sys.stderr, flush=True)
            sys.exit(1)
        return

    # Partition records evenly across n_workers sub-batches
    sub_batches = [records[i::n_workers] for i in range(n_workers)]
    worker_tasks = []
    for w_idx, batch_recs in enumerate(sub_batches):
        if not batch_recs:
            continue
        w_in = Path(f"_gnina_w{w_idx}_in.sdf")
        w_out = Path(f"_gnina_w{w_idx}_out.sdf")
        w_log = Path(f"_gnina_w{w_idx}.log")
        w_in.write_text("".join(batch_recs), encoding="utf-8")
        worker_tasks.append((w_idx, w_in, w_out, w_log))

    results = {}
    try:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(
                    run_single_gnina_worker,
                    w_idx,
                    w_in,
                    w_out,
                    w_log,
                    args,
                    threads_per_worker,
                ): w_idx
                for (w_idx, w_in, w_out, w_log) in worker_tasks
            }
            for fut in as_completed(futures):
                w_idx, w_out, w_log = fut.result()
                results[w_idx] = (w_out, w_log)
    except Exception as e:
        print(str(e), file=sys.stderr, flush=True)
        sys.exit(1)

    # Merge worker outputs in deterministic order
    with (
        open(args.out_sdf, "w", encoding="utf-8") as final_sdf,
        open(args.out_log, "w", encoding="utf-8") as final_log,
    ):
        for w_idx in sorted(results.keys()):
            w_out, w_log = results[w_idx]
            final_sdf.write(w_out.read_text(encoding="utf-8", errors="ignore"))
            final_log.write(
                f"=== WORKER {w_idx} ===\n"
                + w_log.read_text(encoding="utf-8", errors="ignore")
                + "\n"
            )
            w_out.unlink(missing_ok=True)
            w_log.unlink(missing_ok=True)

    for w_idx, w_in, _, _ in worker_tasks:
        w_in.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
