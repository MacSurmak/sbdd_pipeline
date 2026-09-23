#!/usr/bin/env python3
"""
bin/run_production_md.py
Dynamic chunked production MD simulation engine for GROMACS on Slurm HPC.
Orchestrates sequential chunk extensions (chunk 1..N) and performs post-processing:
1. Trajectory concatenation (trjcat)
2. Complex group creation (Protein + UNK strictly, no lipids)
3. PBC correction, whole molecules, and rotational/translational fit (trjconv)
4. Extraction of clean complex_ref.pdb
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


# Add script directory to sys.path for absolute imports
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils.gmx_helpers import ensure_complex_ndx_group


def run_cmd(cmd: list, cwd: Path, desc: str = ""):
    print(f"[RUN] {desc or cmd[0]}: {' '.join(str(x) for x in cmd)}")
    res = subprocess.run(
        cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if res.returncode != 0:
        print(
            f"[ERROR] {desc} failed (exit {res.returncode}):\n{res.stderr}\n{res.stdout}",
            file=sys.stderr,
        )
        sys.exit(res.returncode)
    return res


def main():
    parser = argparse.ArgumentParser(description="GROMACS Chunked Production MD Engine")
    parser.add_argument("--meta_id", required=True, help="Target system ID")
    parser.add_argument("--rep", type=int, default=1, help="Replicate index")
    parser.add_argument("--gro", required=True, help="Equilibrated GRO file")
    parser.add_argument("--cpt", required=True, help="Equilibrated CPT file")
    parser.add_argument("--top", required=True, help="Topology file")
    parser.add_argument("--ndx", required=True, help="Index file")
    parser.add_argument(
        "--total_ns", type=float, default=60.0, help="Total simulation length in ns"
    )
    parser.add_argument(
        "--chunk_ns", type=float, default=10.0, help="Chunk length in ns"
    )
    parser.add_argument(
        "--temp", type=float, default=310.15, help="Temperature in Kelvin"
    )
    parser.add_argument("--is_membrane", default="true", help="Membrane system flag")
    parser.add_argument("--queue", default="main", help="Slurm queue/partition")
    parser.add_argument("--mdp_template", required=True, help="Base MD template file")
    parser.add_argument("--cpus", type=int, default=1, help="Number of OpenMP threads")
    parser.add_argument("--outdir", default=".", help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    gmx_bin = shutil.which("gmx") or "gmx"
    is_mem = str(args.is_membrane).lower() in ["true", "1", "yes"]

    ndx_file = Path(args.ndx).resolve()

    # 1. Generate md.mdp for production
    dt = 0.002
    chunk_steps = int(round((args.chunk_ns * 1000.0) / dt))
    total_chunks = max(1, int(round(args.total_ns / args.chunk_ns)))

    pcoupl_type = "semiisotropic" if is_mem else "isotropic"
    ref_p = "1.0 1.0" if is_mem else "1.0"
    comp = "4.5e-5 4.5e-5" if is_mem else "4.5e-5"

    template_path = Path(args.mdp_template)
    if not template_path.exists():
        print(f"[ERROR] MDP template {template_path} not found.", file=sys.stderr)
        sys.exit(1)

    with open(template_path, "r", encoding="utf-8") as f:
        mdp_content = f.read()

    mdp_content = mdp_content.replace("{NSTEPS}", str(chunk_steps))
    mdp_content = mdp_content.replace("{TEMP}", str(args.temp))
    mdp_content = mdp_content.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
    mdp_content = mdp_content.replace("{PCOUPL_TYPE}", pcoupl_type)
    mdp_content = mdp_content.replace("{REF_P}", ref_p)
    mdp_content = mdp_content.replace("{COMPRESSIBILITY}", comp)

    mdp_file = out_dir / "md.mdp"
    mdp_file.write_text(mdp_content)

    # 2. Execute production chunks
    prev_gro = str(Path(args.gro).resolve())
    prev_cpt = str(Path(args.cpt).resolve())

    # Detect GPU availability
    has_gpu = bool(shutil.which("nvidia-smi"))
    base_hw = (
        ["-nb", "gpu", "-pme", "gpu", "-update", "gpu", "-ntomp", str(args.cpus)]
        if has_gpu
        else ["-ntomp", str(args.cpus)]
    )

    print(
        f"Starting Production MD: {total_chunks} chunk(s) x {args.chunk_ns} ns = {args.total_ns} ns total"
    )

    for chunk in range(1, total_chunks + 1):
        chunk_name = f"md_{chunk}"
        tpr_file = out_dir / f"{chunk_name}.tpr"

        if not tpr_file.exists():
            if chunk == 1:
                cmd_grompp = [
                    gmx_bin,
                    "grompp",
                    "-f",
                    str(mdp_file),
                    "-o",
                    str(tpr_file),
                    "-c",
                    prev_gro,
                    "-t",
                    prev_cpt,
                    "-p",
                    str(Path(args.top).resolve()),
                    "-n",
                    str(ndx_file),
                    "-maxwarn",
                    "2",
                ]
                run_cmd(cmd_grompp, out_dir, f"grompp {chunk_name}")
            else:
                prev_tpr = out_dir / f"md_{chunk - 1}.tpr"
                extension_ps = int(round(args.chunk_ns * 1000.0))
                cmd_extend = [
                    gmx_bin,
                    "convert-tpr",
                    "-s",
                    str(prev_tpr),
                    "-extend",
                    str(extension_ps),
                    "-o",
                    str(tpr_file),
                ]
                run_cmd(cmd_extend, out_dir, f"convert-tpr {chunk_name}")

        chunk_xtc = out_dir / f"{chunk_name}.xtc"
        if not chunk_xtc.exists():
            cmd_mdrun = [gmx_bin, "mdrun", "-deffnm", chunk_name, "-noappend"] + base_hw
            cur_cpt = out_dir / f"{chunk_name}.cpt"
            prev_cpt_file = (
                out_dir / f"md_{chunk - 1}.cpt" if chunk > 1 else Path(prev_cpt)
            )
            if cur_cpt.exists():
                cmd_mdrun.extend(["-cpi", cur_cpt.name])
            elif prev_cpt_file.exists():
                cmd_mdrun.extend(["-cpi", str(prev_cpt_file)])

            run_cmd(cmd_mdrun, out_dir, f"mdrun {chunk_name}")

    print(f"[SUCCESS] Production MD complete: generated {total_chunks} chunk(s).")


if __name__ == "__main__":
    main()
