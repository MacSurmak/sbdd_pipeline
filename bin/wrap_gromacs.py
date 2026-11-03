#!/usr/bin/env python3
import argparse
import sys
import subprocess
import shutil
import re
import socket
import json
import os
from pathlib import Path


def run_mdrun(cmd, task_name):
    """
    Executes GROMACS mdrun.
    Low-level UDP telemetry is handled natively by the containerized /opt/containers/bin/gmx wrapper.
    """
    if "-v" not in cmd:
        cmd.append("-v")

    print(f"Exec: {' '.join(cmd)}")
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        print(f"FATAL: mdrun failed with exit code {res.returncode}", file=sys.stderr)
        sys.exit(res.returncode)


def run_gmx(cmd, input_str=None):
    print(f"Exec: {' '.join(cmd)}")
    res = subprocess.run(cmd, input=input_str, text=True, check=False)
    if res.returncode != 0:
        print(f"FATAL: GROMACS command failed: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(res.returncode)


def cmd_cluster(args):
    tgt = args.xtcs[0]
    if len(args.xtcs) > 1:
        run_gmx(
            ["gmx", "trjcat", "-f"] + args.xtcs + ["-o", "combined_clean.xtc", "-cat"]
        )
        tgt = "combined_clean.xtc"

    run_gmx(["gmx", "make_ndx", "-f", args.ref_pdb, "-o", "index.ndx"], input_str="q\n")
    subprocess.run(
        [
            sys.executable,
            f"{Path(__file__).parent}/gromacs_ndx_fixer.py",
            "--ndx",
            "index.ndx",
            "--action",
            "ensure_complex",
            "--water_names",
            args.water_names,
            "--ion_names",
            args.ion_names,
            "--lipid_names",
            args.lipid_names,
        ],
        check=True,
    )

    grp = args.group
    if grp.lower() == "pocket":
        subprocess.run(
            [
                sys.executable,
                f"{Path(__file__).parent}/gromacs_ndx_fixer.py",
                "--ndx",
                "index.ndx",
                "--pdb",
                args.ref_pdb,
                "--ligand",
                args.ligand,
                "--action",
                "add_pocket",
            ],
            check=True,
        )
        grp = "Pocket"
    elif grp.lower() == "ligand":
        grp = args.ligand

    cmd = [
        "gmx",
        "cluster",
        "-f",
        tgt,
        "-s",
        args.ref_pdb,
        "-n",
        "index.ndx",
        "-method",
        args.method,
        "-cutoff",
        str(args.cutoff),
    ]
    if args.dt and args.dt != "0":
        cmd.extend(["-dt", str(args.dt)])
    if args.max_clusters and int(args.max_clusters) > 0:
        cmd.extend(["-wcl", str(args.max_clusters)])
    cmd.extend(["-cl", "clusters_multi.pdb"])

    run_gmx(cmd, input_str=f"{grp}\nCOMPLEX\n")
    subprocess.run(
        [
            "schrodinger",
            "run",
            "python3",
            f"{Path(__file__).parent}/split_cluster_pdb.py",
            "--input",
            "clusters_multi.pdb",
            "--sys_name",
            args.sys_name,
            "--ligand",
            args.ligand,
            "--max_clusters",
            str(args.max_clusters),
            "--out_receipt",
            args.out_receipt,
        ],
        check=True,
    )


def cmd_trjconv(args):
    tgt = args.xtcs[0]
    if len(args.xtcs) > 1:
        run_gmx(
            ["gmx", "trjcat", "-f"] + args.xtcs + ["-o", "scratch_concat.xtc", "-cat"]
        )
        tgt = "scratch_concat.xtc"

    run_gmx(["gmx", "make_ndx", "-f", args.tpr, "-o", "index.ndx"], input_str="q\n")
    ensure_cmd = [
        sys.executable,
        f"{Path(__file__).parent}/gromacs_ndx_fixer.py",
        "--ndx",
        "index.ndx",
        "--action",
        "ensure_complex",
        "--water_names",
        args.water_names,
        "--ion_names",
        args.ion_names,
        "--lipid_names",
        args.lipid_names,
    ]
    if getattr(args, "dry_top", ""):
        ensure_cmd.extend(["--dry_top", args.dry_top])
    subprocess.run(ensure_cmd, check=True)
    run_gmx(
        [
            "gmx",
            "trjconv",
            "-s",
            args.tpr,
            "-f",
            tgt,
            "-o",
            "pbc_whole.xtc",
            "-pbc",
            "whole",
        ],
        input_str="0\n",
    )
    run_gmx(
        [
            "gmx",
            "trjconv",
            "-s",
            args.tpr,
            "-f",
            "pbc_whole.xtc",
            "-n",
            "index.ndx",
            "-o",
            "pbc_centered.xtc",
            "-center",
            "-pbc",
            "mol",
            "-ur",
            "compact",
        ],
        input_str=f"{args.pbc_center}\n{args.pbc_extract}\n",
    )
    run_gmx(
        [
            "gmx",
            "trjconv",
            "-s",
            args.tpr,
            "-f",
            "pbc_centered.xtc",
            "-n",
            "index.ndx",
            "-o",
            args.out_xtc,
            "-fit",
            "rot+trans",
        ],
        input_str=f"Protein\n{args.pbc_extract}\n",
    )
    run_gmx(
        [
            "gmx",
            "trjconv",
            "-s",
            args.tpr,
            "-f",
            args.out_xtc,
            "-n",
            "index.ndx",
            "-o",
            args.out_pdb,
            "-dump",
            "0",
        ],
        input_str=f"{args.pbc_extract}\n",
    )

    if hasattr(args, "out_ndx") and args.out_ndx:
        run_gmx(
            ["gmx", "make_ndx", "-f", args.out_pdb, "-o", args.out_ndx], input_str="q\n"
        )
        if hasattr(args, "ligand") and args.ligand:
            subprocess.run(
                [
                    sys.executable,
                    f"{Path(__file__).parent}/gromacs_ndx_fixer.py",
                    "--ndx",
                    args.out_ndx,
                    "--pdb",
                    args.out_pdb,
                    "--ligand",
                    args.ligand,
                    "--action",
                    "add_mmpbsa_groups",
                ],
                check=True,
            )

    for tmp in ["pbc_whole.xtc", "pbc_centered.xtc", "scratch_concat.xtc"]:
        if Path(tmp).exists():
            Path(tmp).unlink()


def cmd_equilibration(args):
    if not Path("protocol.txt").exists():
        print("FATAL: protocol.txt playlist not found!", file=sys.stderr)
        sys.exit(1)

    prev_gro = args.input_gro
    prev_cpt = None
    last_stage = None

    with open("protocol.txt", "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split(",")
            if len(parts) != 2:
                print(
                    f"FATAL: Invalid protocol.txt line format (expected type,mdp_path): {line}",
                    file=sys.stderr,
                )
                sys.exit(1)

            p_type, mdp = parts[0].strip(), parts[1].strip()
            stage = Path(mdp).stem
            last_stage = stage

            cmd = [
                "gmx",
                "grompp",
                "-f",
                mdp,
                "-o",
                f"{stage}.tpr",
                "-c",
                prev_gro,
                "-r",
                prev_gro,
                "-p",
                args.topol,
                "-n",
                args.ndx,
                "-maxwarn",
                "0",
            ]

            # Checkpoint file (-t) is passed strictly if it physically exists from a prior dynamical step
            if prev_cpt and Path(prev_cpt).exists() and p_type != "minim":
                cmd.extend(["-t", prev_cpt])

            run_gmx(cmd)

            # Blindly apply exact hardware execution flags from YAML contract
            exec_flags = args.minim_flags if p_type == "minim" else args.eq_flags
            run_cmd = ["gmx", "mdrun", "-deffnm", stage] + exec_flags.split()

            # CPU threading is managed dynamically by Slurm, we append it explicitly
            run_cmd.extend(["-ntomp", str(args.cpus)])

            run_mdrun(run_cmd, f"EQ: {args.sys_name} [{stage}]")

            out_gro = Path(f"{stage}.gro")
            if not out_gro.exists():
                print(
                    f"FATAL: Stage '{stage}' failed to produce coordinate file {out_gro}!",
                    file=sys.stderr,
                )
                sys.exit(1)
            prev_gro = str(out_gro)

            out_cpt = Path(f"{stage}.cpt")
            if out_cpt.exists():
                prev_cpt = str(out_cpt)
            else:
                prev_cpt = None

    if not last_stage:
        print("FATAL: protocol.txt is empty!", file=sys.stderr)
        sys.exit(1)

    shutil.copy2(prev_gro, "eq_final.gro")
    if prev_cpt and Path(prev_cpt).exists():
        shutil.copy2(prev_cpt, "eq_final.cpt")
    else:
        print(
            "FATAL: Multi-stage equilibration finished without generating a velocity checkpoint (.cpt) file! "
            "Production MD cannot be launched without velocities.",
            file=sys.stderr,
        )
        sys.exit(1)

    shutil.copy2(f"{last_stage}.tpr", "eq_final.tpr")


def cmd_run_chunk(args):
    for p in args.posre_files:
        pp = Path(p)
        if pp.is_symlink():
            real = pp.resolve()
            pp.unlink()
            shutil.copy2(real, pp)

    mdp_template = Path("prod.mdp") if Path("prod.mdp").exists() else Path("md.mdp")
    if not mdp_template.exists():
        print(
            "FATAL: Production template (prod.mdp) was not found in staged protocol directory!",
            file=sys.stderr,
        )
        sys.exit(1)

    mdp_content = mdp_template.read_text()
    dt_match = re.search(r"dt\s*=\s*([\d.]+)", mdp_content)
    if not dt_match:
        print(
            "FATAL: Time step 'dt' was not found in production MDP template!",
            file=sys.stderr,
        )
        sys.exit(1)

    dt_ps = float(dt_match.group(1))
    steps = int(round((float(args.chunk_ns) * 1000.0) / dt_ps))
    ext_ps = int(round(float(args.chunk_ns) * 1000.0))

    if int(args.chunk_idx) == 1:
        mdp = re.sub(r"nsteps\s*=.*", f"nsteps = {steps}", mdp_content)
        Path("md_run.mdp").write_text(mdp)

        run_gmx(
            [
                "gmx",
                "grompp",
                "-f",
                "md_run.mdp",
                "-o",
                "out.tpr",
                "-c",
                args.prev_gro,
                "-r",
                args.prev_gro,
                "-t",
                args.prev_cpt,
                "-p",
                args.topol,
                "-n",
                args.ndx,
                "-maxwarn",
                "0",
            ]
        )
        run_cmd = [
            "gmx",
            "mdrun",
            "-deffnm",
            f"md_{args.chunk_idx}",
            "-s",
            "out.tpr",
        ] + args.prod_flags.split()
        run_cmd.extend(["-ntomp", str(args.cpus)])
        run_mdrun(run_cmd, f"MD: {args.sys_name}")
        shutil.copy2("out.tpr", "sys_ref.tpr")
    else:
        run_gmx(
            [
                "gmx",
                "convert-tpr",
                "-s",
                args.prev_tpr,
                "-extend",
                str(ext_ps),
                "-o",
                "out.tpr",
            ]
        )
        run_cmd = [
            "gmx",
            "mdrun",
            "-deffnm",
            f"md_{args.chunk_idx}",
            "-s",
            "out.tpr",
            "-cpi",
            args.prev_cpt,
            "-noappend",
        ] + args.prod_flags.split()
        run_cmd.extend(["-ntomp", str(args.cpus)])
        run_mdrun(run_cmd, f"MD: {args.sys_name}")
        shutil.copy2("out.tpr", "sys_ref.tpr")

    for ext in ["gro", "cpt", "xtc"]:
        for f in Path(".").glob(f"md_{args.chunk_idx}.part*.{ext}"):
            f.rename(f"md_{args.chunk_idx}.{ext}")

    shutil.copy2(f"md_{args.chunk_idx}.gro", "out.gro")
    shutil.copy2(f"md_{args.chunk_idx}.cpt", "out.cpt")


def sanitize_mmpbsa_delta_csv(csv_path: Path):
    """
    Extracts strictly the DELTA binding free energy table from gmx_MMPBSA multi-section CSV,
    eradicating raw protein enthalpies and guaranteeing a single clean table with exact N samples.
    """
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return

    lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    delta_start_idx = None
    for idx, line in enumerate(lines):
        line_upper = line.upper()
        if "DELTA" in line_upper and ("," in line or "COMPLEX" in line_upper):
            delta_start_idx = idx
            break

    if delta_start_idx is None:
        for idx in range(len(lines) - 1, -1, -1):
            if "Frame" in lines[idx] and (
                "TOTAL" in lines[idx] or "DELTA" in lines[idx]
            ):
                delta_start_idx = idx
                break

    if delta_start_idx is not None:
        if "Frame" not in lines[delta_start_idx]:
            for i in range(delta_start_idx, min(delta_start_idx + 5, len(lines))):
                if "Frame" in lines[i]:
                    delta_start_idx = i
                    break

        clean_lines = []
        for line in lines[delta_start_idx:]:
            l_clean = line.strip()
            if not l_clean or l_clean.startswith("#"):
                if clean_lines and len(clean_lines) > 5:
                    break
                continue
            clean_lines.append(l_clean)

        if clean_lines:
            raw_copy = csv_path.with_name("raw_" + csv_path.name)
            shutil.copy2(csv_path, raw_copy)
            csv_path.write_text("\n".join(clean_lines) + "\n", encoding="utf-8")
            print(
                f"[MMPBSA] Sanitized DELTA section in {csv_path.name} ({len(clean_lines)-1} frames retained).",
                flush=True,
            )


def cmd_mmpbsa(args):
    import MDAnalysis as mda

    u = mda.Universe(args.ref_pdb, args.clean_xtc)
    n_frames = len(u.trajectory)
    if n_frames < 2:
        raise ValueError(
            f"FATAL: Trajectory {args.clean_xtc} has fewer than 2 frames ({n_frames}). Cannot evaluate autocorrelation or sampling stride!"
        )

    dt_ps = float(u.trajectory[1].time - u.trajectory[0].time)
    if dt_ps <= 0:
        raise ValueError(
            f"FATAL: Non-physical trajectory timestamp interval detected: dt = {dt_ps} ps in {args.clean_xtc}"
        )

    dt_ns = dt_ps / 1000.0
    sample_interval_ps = float(args.sample_interval_ps)
    if sample_interval_ps < dt_ps:
        raise ValueError(
            f"FATAL: Requested MM-PBSA sample_interval_ps ({sample_interval_ps:.2f} ps) is smaller than trajectory "
            f"frame recording interval ({dt_ps:.2f} ps)! Cannot sample faster than trajectory frequency."
        )

    step_frames = max(1, int(round(sample_interval_ps / dt_ps)))
    start_frame = min(n_frames, max(1, int(round(float(args.start_ns) / dt_ns)) + 1))

    actual_interval_ps = step_frames * dt_ps
    print(
        f"[MMPBSA] Sampling setup: dt_traj={dt_ps:.2f} ps, target_interval={sample_interval_ps:.2f} ps, "
        f"step_frames={step_frames} (actual_interval={actual_interval_ps:.2f} ps), start_frame={start_frame}/{n_frames}",
        flush=True,
    )

    # 1. Resolve calculation template declaratively from assets/mmpbsa/
    tmpl_dir = Path(args.templates_dir).resolve()
    if not tmpl_dir.exists():
        raise FileNotFoundError(
            f"FATAL: MM-PBSA templates directory not found: {tmpl_dir}"
        )

    method_req = str(args.method).strip().lower()
    if method_req == "auto":
        target_template_name = "pb_mem.in" if args.membrane else "gb.in"
    elif method_req == "pb":
        target_template_name = "pb_mem.in" if args.membrane else "pb_solv.in"
    elif method_req == "gb":
        target_template_name = "gb.in"
    else:
        raise ValueError(f"FATAL: Unsupported MM-PBSA method '{method_req}'!")

    primary_tmpl_path = tmpl_dir / target_template_name
    if not primary_tmpl_path.exists():
        raise FileNotFoundError(
            f"FATAL: Archetype template {target_template_name} not found in {tmpl_dir}!"
        )

    # Format primary calculation input file
    tmpl_text = primary_tmpl_path.read_text(encoding="utf-8")
    formatted_in = (
        tmpl_text.replace("{SYS_NAME}", "Protein_Ligand")
        .replace("{START_FRAME}", str(start_frame))
        .replace("{STEP_FRAMES}", str(step_frames))
        .replace("{TEMP}", str(args.temp))
        .replace("{SALT_CONC}", str(args.salt_conc))
        .replace("{IGB}", str(args.gb_model))
        .replace("{MEM_THICKNESS}", str(args.membrane_thickness))
        .replace("{EMEM}", str(args.emem))
        .replace("{INDI}", str(args.indi))
        .replace("{FILL_RATIO}", str(args.fillratio))
        .replace("{IPB}", "1")
        .replace("{INP}", "2")
    )
    Path("mmpbsa_run.in").write_text(formatted_in, encoding="utf-8")

    mpirun = (
        ["mpirun", "--oversubscribe", "-np", str(args.cpus), "gmx_MMPBSA", "MPI"]
        if int(args.cpus) > 1
        else ["gmx_MMPBSA"]
    )

    run_gmx(
        mpirun
        + [
            "-O",
            "-i",
            "mmpbsa_run.in",
            "-cs",
            args.ref_pdb,
            "-ci",
            args.ndx,
            "-cg",
            "MMPBSA_REC",
            "MMPBSA_LIG",
            "-ct",
            args.clean_xtc,
            "-cp",
            args.topol,
            "-nogui",
            "-o",
            "thermodynamics_mmpbsa_pb.dat",
            "-eo",
            "thermodynamics_mmpbsa_pb.csv",
        ]
    )

    # Sanitize CSV: extract strictly the DELTA table so downstream tools receive clean binding energies
    sanitize_mmpbsa_delta_csv(Path("thermodynamics_mmpbsa_pb.csv"))

    # 2. Per-residue decomposition (strictly GB-based)
    if args.decomp:
        decomp_tmpl_path = tmpl_dir / "decomp_gb.in"
        if not decomp_tmpl_path.exists():
            raise FileNotFoundError(
                f"FATAL: Decomposition template decomp_gb.in not found in {tmpl_dir}!"
            )

        decomp_text = decomp_tmpl_path.read_text(encoding="utf-8")
        formatted_decomp = (
            decomp_text.replace("{SYS_NAME}", "Protein_Ligand")
            .replace("{START_FRAME}", str(start_frame))
            .replace("{STEP_FRAMES}", str(step_frames))
            .replace("{TEMP}", str(args.temp))
            .replace("{SALT_CONC}", str(args.salt_conc))
            .replace("{IGB}", str(args.gb_model))
            .replace("{IDECOMP}", "2")
        )
        Path("mmpbsa_decomp.in").write_text(formatted_decomp, encoding="utf-8")

        run_gmx(
            mpirun
            + [
                "-O",
                "-prefix",
                "_GMXMMPBSA_GB_",
                "-i",
                "mmpbsa_decomp.in",
                "-cs",
                args.ref_pdb,
                "-ci",
                args.ndx,
                "-cg",
                "MMPBSA_REC",
                "MMPBSA_LIG",
                "-ct",
                args.clean_xtc,
                "-cp",
                args.topol,
                "-nogui",
                "-o",
                "thermodynamics_mmpbsa_gb.dat",
                "-deo",
                "decomposition_mmpbsa_gb.csv",
            ]
        )


def cmd_trjcat(args):
    if len(args.xtcs) > 1:
        run_gmx(["gmx", "trjcat", "-f"] + args.xtcs + ["-cat", "-o", args.out_xtc])
    else:
        shutil.copy2(args.xtcs[0], args.out_xtc)


def cmd_ndx_posre(args):
    run_gmx(
        ["gmx", "make_ndx", "-f", args.input_gro, "-o", "index.ndx"], input_str="q\n"
    )
    subprocess.run(
        [
            sys.executable,
            f"{Path(__file__).parent}/gromacs_ndx_fixer.py",
            "--ndx",
            "index.ndx",
            "--action",
            "add_groups",
            "--water_names",
            args.water_names,
            "--ion_names",
            args.ion_names,
        ],
        check=True,
    )
    cmd_gen = [
        sys.executable,
        f"{Path(__file__).parent}/gromacs_ndx_fixer.py",
        "--ndx",
        "index.ndx",
        "--top",
        args.topol,
        "--forces",
        args.forces,
        "--action",
        "generate_posre",
        "--water_names",
        args.water_names,
        "--ion_names",
        args.ion_names,
    ]
    if hasattr(args, "prod_fc") and args.prod_fc and int(args.prod_fc) > 0:
        cmd_gen.extend(
            [
                "--pdb",
                args.input_gro,
                "--prod_fc",
                str(args.prod_fc),
                "--prod_group",
                getattr(args, "prod_group", "C-alpha"),
            ]
        )
    subprocess.run(cmd_gen, check=True)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    p_cluster = subparsers.add_parser("cluster")
    p_cluster.add_argument("--xtcs", nargs="+", required=True)
    p_cluster.add_argument("--ref_pdb", required=True)
    p_cluster.add_argument("--cutoff", required=True)
    p_cluster.add_argument("--method", required=True)
    p_cluster.add_argument("--group", required=True)
    p_cluster.add_argument("--ligand", required=True)
    p_cluster.add_argument("--max_clusters", required=True)
    p_cluster.add_argument("--sys_name", required=True)
    p_cluster.add_argument("--dt")
    p_cluster.add_argument("--water_names", required=True)
    p_cluster.add_argument("--ion_names", required=True)
    p_cluster.add_argument("--lipid_names", required=True)
    p_cluster.add_argument("--out_receipt", default="cluster_receipt.json")

    p_trjconv = subparsers.add_parser("trjconv")
    p_trjconv.add_argument("--xtcs", nargs="+", required=True)
    p_trjconv.add_argument("--tpr", required=True)
    p_trjconv.add_argument("--pbc_center", required=True)
    p_trjconv.add_argument("--pbc_extract", required=True)
    p_trjconv.add_argument("--out_xtc", required=True)
    p_trjconv.add_argument("--out_pdb", required=True)
    p_trjconv.add_argument("--out_ndx", default="")
    p_trjconv.add_argument("--ligand", default="")
    p_trjconv.add_argument("--dry_top", default="")
    p_trjconv.add_argument("--water_names", required=True)
    p_trjconv.add_argument("--ion_names", required=True)
    p_trjconv.add_argument("--lipid_names", required=True)

    p_eq = subparsers.add_parser("equilibration")
    p_eq.add_argument("--input_gro", required=True)
    p_eq.add_argument("--topol", required=True)
    p_eq.add_argument("--ndx", required=True)
    p_eq.add_argument("--cpus", required=True)
    p_eq.add_argument("--posre_files", nargs="+", required=True)
    p_eq.add_argument("--sys_name", required=True)
    p_eq.add_argument("--minim_flags", required=True)
    p_eq.add_argument("--eq_flags", required=True)

    p_chunk = subparsers.add_parser("run_chunk")
    p_chunk.add_argument("--prev_gro", required=True)
    p_chunk.add_argument("--prev_cpt", required=True)
    p_chunk.add_argument("--prev_tpr", required=True)
    p_chunk.add_argument("--topol", required=True)
    p_chunk.add_argument("--ndx", required=True)
    p_chunk.add_argument("--chunk_ns", required=True)
    p_chunk.add_argument("--chunk_idx", required=True)
    p_chunk.add_argument("--cpus", required=True)
    p_chunk.add_argument("--posre_files", nargs="+", required=True)
    p_chunk.add_argument("--sys_name", required=True)
    p_chunk.add_argument("--prod_flags", required=True)

    p_mmpbsa = subparsers.add_parser("mmpbsa")
    p_mmpbsa.add_argument("--clean_xtc", required=True)
    p_mmpbsa.add_argument("--ref_pdb", required=True)
    p_mmpbsa.add_argument("--topol", required=True)
    p_mmpbsa.add_argument("--ndx", required=True)
    p_mmpbsa.add_argument("--ligand", required=True)
    p_mmpbsa.add_argument("--templates_dir", required=True)
    p_mmpbsa.add_argument("--start_ns", required=True)
    p_mmpbsa.add_argument("--sample_interval_ps", required=True)
    p_mmpbsa.add_argument("--temp", required=True)
    p_mmpbsa.add_argument("--salt_conc", required=True)
    p_mmpbsa.add_argument("--cpus", required=True)
    p_mmpbsa.add_argument("--method", default="auto")
    p_mmpbsa.add_argument("--gb_model", default="2")
    p_mmpbsa.add_argument("--emem", default="7.0")
    p_mmpbsa.add_argument("--indi", default="4.0")
    p_mmpbsa.add_argument("--fillratio", default="1.25")
    p_mmpbsa.add_argument("--membrane", action="store_true")
    p_mmpbsa.add_argument("--membrane_thickness", default="28.0")
    p_mmpbsa.add_argument("--decomp", action="store_true")

    p_cat = subparsers.add_parser("trjcat")
    p_cat.add_argument("--xtcs", nargs="+", required=True)
    p_cat.add_argument("--out_xtc", required=True)

    p_ndx = subparsers.add_parser("ndx_posre")
    p_ndx.add_argument("--input_gro", required=True)
    p_ndx.add_argument("--topol", required=True)
    p_ndx.add_argument("--forces", required=True)
    p_ndx.add_argument("--water_names", required=True)
    p_ndx.add_argument("--ion_names", required=True)
    p_ndx.add_argument("--lipid_names", required=True)
    p_ndx.add_argument("--prod_fc", type=int, default=0)
    p_ndx.add_argument("--prod_group", default="C-alpha")

    args = parser.parse_args()
    if args.subcommand == "cluster":
        cmd_cluster(args)
    elif args.subcommand == "trjconv":
        cmd_trjconv(args)
    elif args.subcommand == "equilibration":
        cmd_equilibration(args)
    elif args.subcommand == "run_chunk":
        cmd_run_chunk(args)
    elif args.subcommand == "mmpbsa":
        cmd_mmpbsa(args)
    elif args.subcommand == "trjcat":
        cmd_trjcat(args)
    elif args.subcommand == "ndx_posre":
        cmd_ndx_posre(args)


if __name__ == "__main__":
    main()
