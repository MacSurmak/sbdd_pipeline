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
    parser = argparse.ArgumentParser(description="Generate GROMACS MDP files from templates")
    parser.add_argument("--templates_dir", required=True, help="Directory containing template MDPs")
    parser.add_argument("--outdir", required=True, help="Output directory for generated MDPs")
    parser.add_argument("--is_membrane", default="false", help="Whether system is a membrane system")
    parser.add_argument("--temperature", type=float, default=310.15, help="Simulation temperature in K")
    parser.add_argument("--nvt_short_ps", type=float, default=50.0, help="NVT short duration in ps")
    parser.add_argument("--nvt_long_ps", type=float, default=500.0, help="NVT long duration in ps")
    parser.add_argument("--npt_berendsen_ps", type=float, default=1000.0, help="NPT Berendsen duration in ps")
    parser.add_argument("--npt_pr_ps", type=float, default=1000.0, help="NPT PR duration in ps")
    args = parser.parse_args()

    tpl_dir = Path(args.templates_dir)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    is_mem = str(args.is_membrane).lower() in ["true", "1", "yes"]
    temp = args.temperature

    pcoupl_type = "semiisotropic" if is_mem else "isotropic"
    ref_p = "1.0 1.0" if is_mem else "1.0"
    compressibility = "4.5e-5 4.5e-5" if is_mem else "4.5e-5"

    dt_minim = 0.01
    dt_nvt = 0.001
    dt_npt = 0.002
    dt_md = 0.002

    # 1. Minimization
    minim_tpl = (tpl_dir / "minim.mdp").read_text()
    minim_content = minim_tpl.replace("{DEFINE}", "define = -DPOSRES -DPOSRES_LIG")
    minim_content = minim_content.replace("{EMSTEP}", str(dt_minim))
    minim_content = minim_content.replace("{CONSTRAINTS}", "h-bonds")
    (out_dir / "minim.mdp").write_text(minim_content)

    # 2. NVT Stages
    nvt_tpl = (tpl_dir / "nvt.mdp").read_text()

    def build_nvt(name, duration_ps, define_str, cont, gen_vel):
        steps = int(round(duration_ps / dt_nvt))
        c = nvt_tpl.replace("{DEFINE}", define_str)
        c = c.replace("{DT}", str(dt_nvt))
        c = c.replace("{NSTEPS}", str(steps))
        c = c.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
        c = c.replace("{TEMP}", str(temp))
        c = c.replace("{CONTINUATION}", cont)
        c = c.replace("{GEN_VEL}", gen_vel)
        (out_dir / f"{name}.mdp").write_text(c)

    build_nvt("nvt_short", args.nvt_short_ps, "define = -DPOSRES -DPOSRES_LIG", "no", "yes")
    build_nvt("nvt_long", args.nvt_long_ps, "define = -DPOSRES_250 -DPOSRES_LIG_250", "yes", "no")
    build_nvt("nvt", args.nvt_short_ps, "define = -DPOSRES -DPOSRES_LIG", "no", "yes")

    # 3. NPT Stages
    npt_tpl = (tpl_dir / "npt.mdp").read_text()

    def build_npt(name, duration_ps, define_str, pcoupl, tau_p):
        steps = int(round(duration_ps / dt_npt))
        c = npt_tpl.replace("{DEFINE}", define_str)
        c = c.replace("{NSTEPS}", str(steps))
        c = c.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
        c = c.replace("{TEMP}", str(temp))
        c = c.replace("{PCOUPL}", pcoupl)
        c = c.replace("{PCOUPL_TYPE}", pcoupl_type)
        c = c.replace("{TAU_P}", str(tau_p))
        c = c.replace("{REF_P}", ref_p)
        c = c.replace("{COMPRESSIBILITY}", compressibility)
        (out_dir / f"{name}.mdp").write_text(c)

    build_npt("npt_berendsen", args.npt_berendsen_ps, "define = -DPOSRES_250 -DPOSRES_LIG_250", "berendsen", 1.0)
    build_npt("npt_pr", args.npt_pr_ps, "define = -DPOSRES_50 -DPOSRES_LIG_50", "Parrinello-Rahman", 5.0)
    build_npt("npt", args.npt_berendsen_ps, "define = -DPOSRES_250 -DPOSRES_LIG_250", "Parrinello-Rahman", 5.0)

    # 4. Production MD
    md_tpl = (tpl_dir / "md.mdp").read_text()
    md_content = md_tpl.replace("{NSTEPS}", "5000000")
    md_content = md_content.replace("{NSTXOUT}", "5000")
    md_content = md_content.replace("{NSTENERGY}", "5000")
    md_content = md_content.replace("{TC_GRPS}", "PROT_LIG SOL_ION")
    md_content = md_content.replace("{TEMP}", str(temp))
    md_content = md_content.replace("{PCOUPL_TYPE}", pcoupl_type)
    md_content = md_content.replace("{REF_P}", ref_p)
    md_content = md_content.replace("{COMPRESSIBILITY}", compressibility)
    (out_dir / "md.mdp").write_text(md_content)

    # 5. Playlists
    (out_dir / "stages_list.txt").write_text("minim.mdp\nnvt_short.mdp\nnvt_long.mdp\nnpt_berendsen.mdp\nnpt_pr.mdp\n")
    (out_dir / "protocol.txt").write_text("minim.mdp\nnvt.mdp\nnpt.mdp\n")

    print(f"Generated complete MDP suite in {out_dir}")


if __name__ == "__main__":
    main()
