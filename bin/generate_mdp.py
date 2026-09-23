#!/usr/bin/env python3
"""
bin/generate_mdp.py
Generates GROMACS MDP files from base templates using dynamic parameters.
"""
import argparse
from pathlib import Path


def ps_to_steps(ps, dt):
    return str(int(round(float(ps) / float(dt))))


def main():
    parser = argparse.ArgumentParser(description="Generate parameterized MDP files")
    parser.add_argument("--templates_dir", required=True)
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--is_membrane", required=True)
    parser.add_argument("--temperature", required=True)
    parser.add_argument("--nvt_short_ps", required=True)
    parser.add_argument("--nvt_long_ps", required=True)
    parser.add_argument("--npt_berendsen_ps", required=True)
    parser.add_argument("--npt_pr_ps", required=True)
    args = parser.parse_args()

    is_mem = str(args.is_membrane).lower() in ["true", "1", "yes"]
    temp = args.temperature

    nvt_short_steps = ps_to_steps(args.nvt_short_ps, 0.001)
    nvt_long_steps = ps_to_steps(args.nvt_long_ps, 0.002)
    npt_berendsen_steps = ps_to_steps(args.npt_berendsen_ps, 0.002)
    npt_pr_steps = ps_to_steps(args.npt_pr_ps, 0.002)

    pcoupl_type = "semiisotropic" if is_mem else "isotropic"
    ref_p = "1.0 1.0" if is_mem else "1.0"
    comp = "4.5e-5 4.5e-5" if is_mem else "4.5e-5"

    base_replace = {
        "{TEMP}": temp,
        "{TC_GRPS}": "PROT_LIG SOL_ION",
        "{PCOUPL_TYPE}": pcoupl_type,
        "{REF_P}": ref_p,
        "{COMPRESSIBILITY}": comp,
    }

    if is_mem:
        stages = [
            (
                "emA",
                "minim.mdp",
                {
                    "{DEFINE}": "-DPOSRES -DPOSRES_LIG",
                    "{EMSTEP}": "0.005",
                    "{CONSTRAINTS}": "h-bonds",
                },
            ),
            (
                "emB",
                "minim.mdp",
                {"{DEFINE}": "", "{EMSTEP}": "0.001", "{CONSTRAINTS}": "h-bonds"},
            ),
            (
                "emC",
                "minim.mdp",
                {"{DEFINE}": "", "{EMSTEP}": "0.0001", "{CONSTRAINTS}": "none"},
            ),
            (
                "nvt_short",
                "nvt.mdp",
                {
                    "{DEFINE}": "-DPOSRES -DPOSRES_LIG",
                    "{NSTEPS}": nvt_short_steps,
                    "{DT}": "0.001",
                    "{GEN_VEL}": "yes",
                    "{CONTINUATION}": "no",
                },
            ),
            (
                "nvt_long",
                "nvt.mdp",
                {
                    "{DEFINE}": "-DPOSRES -DPOSRES_LIG",
                    "{NSTEPS}": nvt_long_steps,
                    "{DT}": "0.002",
                    "{GEN_VEL}": "no",
                    "{CONTINUATION}": "yes",
                },
            ),
            (
                "npt_berendsen",
                "npt.mdp",
                {
                    "{DEFINE}": "-DPOSRES -DPOSRES_LIG",
                    "{NSTEPS}": npt_berendsen_steps,
                    "{PCOUPL}": "C-rescale",
                    "{TAU_P}": "5.0",
                },
            ),
            (
                "npt_pr_250",
                "npt.mdp",
                {
                    "{DEFINE}": "-DPOSRES_250 -DPOSRES_LIG_250",
                    "{NSTEPS}": npt_pr_steps,
                    "{PCOUPL}": "Parrinello-Rahman",
                    "{TAU_P}": "5.0",
                },
            ),
            (
                "npt_pr_50",
                "npt.mdp",
                {
                    "{DEFINE}": "-DPOSRES_50 -DPOSRES_LIG_50",
                    "{NSTEPS}": npt_pr_steps,
                    "{PCOUPL}": "Parrinello-Rahman",
                    "{TAU_P}": "5.0",
                },
            ),
            (
                "npt_unrestrained",
                "npt.mdp",
                {
                    "{DEFINE}": "",
                    "{NSTEPS}": npt_pr_steps,
                    "{PCOUPL}": "Parrinello-Rahman",
                    "{TAU_P}": "5.0",
                },
            ),
        ]
    else:
        stages = [
            (
                "emA",
                "minim.mdp",
                {
                    "{DEFINE}": "-DPOSRES -DPOSRES_LIG",
                    "{EMSTEP}": "0.01",
                    "{CONSTRAINTS}": "h-bonds",
                },
            ),
            (
                "nvt_short",
                "nvt.mdp",
                {
                    "{DEFINE}": "-DPOSRES -DPOSRES_LIG",
                    "{NSTEPS}": nvt_short_steps,
                    "{DT}": "0.001",
                    "{GEN_VEL}": "yes",
                    "{CONTINUATION}": "no",
                },
            ),
            (
                "npt_pr_1000",
                "npt.mdp",
                {
                    "{DEFINE}": "-DPOSRES -DPOSRES_LIG",
                    "{NSTEPS}": npt_pr_steps,
                    "{PCOUPL}": "Parrinello-Rahman",
                    "{TAU_P}": "5.0",
                },
            ),
            (
                "npt_unrestrained",
                "npt.mdp",
                {
                    "{DEFINE}": "",
                    "{NSTEPS}": npt_pr_steps,
                    "{PCOUPL}": "Parrinello-Rahman",
                    "{TAU_P}": "5.0",
                },
            ),
        ]

    out_dir = Path(args.outdir)
    tpl_dir = Path(args.templates_dir)

    for name, tpl, macros in stages:
        with open(tpl_dir / tpl) as f:
            c = f.read()
        for k, v in base_replace.items():
            c = c.replace(k, v)
        for k, v in macros.items():
            c = c.replace(k, v)
        with open(out_dir / f"{name}.mdp", "w") as out:
            out.write(c)

    with open(out_dir / "stages_list.txt", "w") as f:
        for name, _, _ in stages:
            f.write(f"{name}\n")


if __name__ == "__main__":
    main()
