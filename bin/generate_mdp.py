#!/usr/bin/env python3
"""
bin/generate_mdp.py
Generates customized GROMACS MDP files from templates in assets/mdp.
Fills placeholders for temperature, steps, coupling, and restraints without silent fallbacks.
"""
import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Generate pure, system-specific GROMACS MDP protocol suite"
    )
    parser.add_argument(
        "--templates_dir", required=True, help="Directory containing template MDPs"
    )
    parser.add_argument(
        "--outdir", required=True, help="Output directory for generated MDPs"
    )
    parser.add_argument(
        "--is_membrane",
        required=True,
        help="Whether system is a membrane system (true/false)",
    )
    parser.add_argument(
        "--temperature", type=float, required=True, help="Simulation temperature in K"
    )
    parser.add_argument(
        "--chunk_ns",
        type=float,
        required=True,
        help="Production chunk size in ns directly from md.yaml",
    )
    parser.add_argument("--dt_ps", required=True)
    parser.add_argument("--compressibility", required=True)
    parser.add_argument("--tau_p", required=True)
    parser.add_argument("--disp_corr", required=True)
    parser.add_argument("--lincs_order", type=int, required=True)
    parser.add_argument("--lincs_iter", type=int, required=True)
    parser.add_argument(
        "--prod_restraints_fc",
        default="",
        help="Optional force constant for production restraints (e.g. 50)",
    )
    args = parser.parse_args()

    tpl_dir = Path(args.templates_dir)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Purge any old, stale MDP files and lists in the target directory
    for stale_file in out_dir.glob("*.mdp"):
        stale_file.unlink()
    for stale_list in ["protocol.txt", "stages_list.txt"]:
        (out_dir / stale_list).unlink(missing_ok=True)

    is_mem = str(args.is_membrane).lower() in ["true", "1", "yes"]
    temp = args.temperature
    dt_prod = float(args.dt_ps)
    lincs_order_str = str(args.lincs_order)
    lincs_iter_str = str(args.lincs_iter)

    pcoupl_type = "semiisotropic" if is_mem else "isotropic"
    ref_p = "1.0 1.0" if is_mem else "1.0"

    comp_base = args.compressibility
    compressibility = f"{comp_base} {comp_base}" if is_mem else comp_base

    disp_corr = str(args.disp_corr).strip()

    em_tpl = (tpl_dir / "em.mdp").read_text()
    nvt_tpl = (tpl_dir / "nvt.mdp").read_text()
    npt_tpl = (tpl_dir / "npt.mdp").read_text()
    prod_tpl = (tpl_dir / "prod.mdp").read_text()

    # 1. EM: Energy Minimization (Steepest descent with strong restraints)
    em_content = em_tpl.replace("{DEFINE}", "-DPOSRES -DPOSRES_LIG")
    em_content = em_content.replace("{DISPCORR}", disp_corr)
    (out_dir / "em.mdp").write_text(em_content)

    # 2. NVT: Thermalization / Heating
    nvt_content = nvt_tpl.replace("{DEFINE}", "-DPOSRES -DPOSRES_LIG")
    nvt_content = nvt_content.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
    nvt_content = nvt_content.replace("{TEMP}", str(temp))
    nvt_content = nvt_content.replace("{DISPCORR}", disp_corr)
    nvt_content = nvt_content.replace("{LINCS_ORDER}", lincs_order_str)
    nvt_content = nvt_content.replace("{LINCS_ITER}", lincs_iter_str)
    (out_dir / "nvt.mdp").write_text(nvt_content)

    # 3. NPT Stages & Playlist
    if is_mem:
        # Stepwise Membrane Relaxation:
        # Stage 1: Stochastic cell rescaling (C-rescale) relaxation to pack lipid area without oscillations
        npt_relax = npt_tpl.replace("{DEFINE}", "-DPOSRES_250 -DPOSRES_LIG_250")
        npt_relax = npt_relax.replace("{DT}", "0.002")
        npt_relax = npt_relax.replace("{NSTEPS}", "500000")  # 1 ns at 2 fs
        npt_relax = npt_relax.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
        npt_relax = npt_relax.replace("{TEMP}", str(temp))
        npt_relax = npt_relax.replace("{PCOUPL}", "C-rescale")
        npt_relax = npt_relax.replace("{PCOUPL_TYPE}", pcoupl_type)
        npt_relax = npt_relax.replace("{TAU_P}", "2.0")
        npt_relax = npt_relax.replace("{REF_P}", ref_p)
        npt_relax = npt_relax.replace("{COMPRESSIBILITY}", compressibility)
        npt_relax = npt_relax.replace("{DISPCORR}", disp_corr)
        npt_relax = npt_relax.replace("{LINCS_ORDER}", lincs_order_str)
        npt_relax = npt_relax.replace("{LINCS_ITER}", lincs_iter_str)
        (out_dir / "npt_relax.mdp").write_text(npt_relax)

        # Stage 2: Parrinello-Rahman density equilibration at target production dt_ps
        npt_equil_steps = int(round(1000.0 / dt_prod))  # 1 ns
        npt_equil = npt_tpl.replace("{DEFINE}", "-DPOSRES_50 -DPOSRES_LIG_50")
        npt_equil = npt_equil.replace("{DT}", str(args.dt_ps))
        npt_equil = npt_equil.replace("{NSTEPS}", str(npt_equil_steps))
        npt_equil = npt_equil.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
        npt_equil = npt_equil.replace("{TEMP}", str(temp))
        npt_equil = npt_equil.replace("{PCOUPL}", "Parrinello-Rahman")
        npt_equil = npt_equil.replace("{PCOUPL_TYPE}", pcoupl_type)
        npt_equil = npt_equil.replace("{TAU_P}", str(args.tau_p))
        npt_equil = npt_equil.replace("{REF_P}", ref_p)
        npt_equil = npt_equil.replace("{COMPRESSIBILITY}", compressibility)
        npt_equil = npt_equil.replace("{DISPCORR}", disp_corr)
        npt_equil = npt_equil.replace("{LINCS_ORDER}", lincs_order_str)
        npt_equil = npt_equil.replace("{LINCS_ITER}", lincs_iter_str)
        (out_dir / "npt_equil.mdp").write_text(npt_equil)

        playlist = [
            "minim,em.mdp",
            "eq,nvt.mdp",
            "eq,npt_relax.mdp",
            "eq,npt_equil.mdp",
        ]
    else:
        # Water/Solvent Relaxation: 2-stage density equilibration to eliminate pressure shock
        # Stage 1: C-rescale relaxation to gently compress solvent to 1.0 g/cm3 without box oscillation
        npt_relax = npt_tpl.replace("{DEFINE}", "-DPOSRES_250 -DPOSRES_LIG_250")
        npt_relax = npt_relax.replace("{DT}", "0.002")
        npt_relax = npt_relax.replace("{NSTEPS}", "100000")  # 200 ps at 2 fs
        npt_relax = npt_relax.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
        npt_relax = npt_relax.replace("{TEMP}", str(temp))
        npt_relax = npt_relax.replace("{PCOUPL}", "C-rescale")
        npt_relax = npt_relax.replace("{PCOUPL_TYPE}", pcoupl_type)
        npt_relax = npt_relax.replace("{TAU_P}", "2.0")
        npt_relax = npt_relax.replace("{REF_P}", ref_p)
        npt_relax = npt_relax.replace("{COMPRESSIBILITY}", compressibility)
        npt_relax = npt_relax.replace("{DISPCORR}", disp_corr)
        npt_relax = npt_relax.replace("{LINCS_ORDER}", lincs_order_str)
        npt_relax = npt_relax.replace("{LINCS_ITER}", lincs_iter_str)
        (out_dir / "npt_relax.mdp").write_text(npt_relax)

        # Stage 2: Parrinello-Rahman true isothermal-isobaric ensemble equilibration at target production dt_ps
        npt_equil_steps = int(round(800.0 / dt_prod))  # 800 ps
        npt_equil = npt_tpl.replace("{DEFINE}", "-DPOSRES_50 -DPOSRES_LIG_50")
        npt_equil = npt_equil.replace("{DT}", str(args.dt_ps))
        npt_equil = npt_equil.replace("{NSTEPS}", str(npt_equil_steps))
        npt_equil = npt_equil.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
        npt_equil = npt_equil.replace("{TEMP}", str(temp))
        npt_equil = npt_equil.replace("{PCOUPL}", "Parrinello-Rahman")
        npt_equil = npt_equil.replace("{PCOUPL_TYPE}", pcoupl_type)
        npt_equil = npt_equil.replace("{TAU_P}", str(args.tau_p))
        npt_equil = npt_equil.replace("{REF_P}", ref_p)
        npt_equil = npt_equil.replace("{COMPRESSIBILITY}", compressibility)
        npt_equil = npt_equil.replace("{DISPCORR}", disp_corr)
        npt_equil = npt_equil.replace("{LINCS_ORDER}", lincs_order_str)
        npt_equil = npt_equil.replace("{LINCS_ITER}", lincs_iter_str)
        (out_dir / "npt_equil.mdp").write_text(npt_equil)

        playlist = [
            "minim,em.mdp",
            "eq,nvt.mdp",
            "eq,npt_relax.mdp",
            "eq,npt_equil.mdp",
        ]

    # 4. Production Template (Directly synchronized with chunk_size_ns and dt_ps from md.yaml)
    prod_steps = int(round((args.chunk_ns * 1000.0) / dt_prod))
    nstxout = max(1, min(int(round(10.0 / dt_prod)), prod_steps // 10))
    nstenergy = max(1, min(int(round(2.0 / dt_prod)), prod_steps // 10))

    prod_define = ""
    if args.prod_restraints_fc and float(args.prod_restraints_fc) > 0:
        prod_define = f"define                   = -DPOSRES_PROD_{int(float(args.prod_restraints_fc))}"

    prod_content = prod_tpl.replace("{DEFINE}", prod_define)
    prod_content = prod_content.replace("{DT}", str(args.dt_ps))
    prod_content = prod_content.replace("{NSTEPS}", str(prod_steps))
    prod_content = prod_content.replace("{NSTXOUT}", str(nstxout))
    prod_content = prod_content.replace("{NSTENERGY}", str(nstenergy))
    prod_content = prod_content.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
    prod_content = prod_content.replace("{TEMP}", str(temp))
    prod_content = prod_content.replace("{PCOUPL_TYPE}", pcoupl_type)
    prod_content = prod_content.replace("{TAU_P}", str(args.tau_p))
    prod_content = prod_content.replace("{REF_P}", ref_p)
    prod_content = prod_content.replace("{COMPRESSIBILITY}", compressibility)
    prod_content = prod_content.replace("{DISPCORR}", disp_corr)
    prod_content = prod_content.replace("{LINCS_ORDER}", lincs_order_str)
    prod_content = prod_content.replace("{LINCS_ITER}", lincs_iter_str)
    (out_dir / "prod.mdp").write_text(prod_content)

    # 5. Master Equilibration Playlist
    (out_dir / "protocol.txt").write_text("\n".join(playlist) + "\n")

    print(
        f"Generated clean {('membrane' if is_mem else 'water')} MD protocol suite ({len(playlist) + 1} files) in {out_dir}"
    )


if __name__ == "__main__":
    main()
