/*
 * modules/gromacs/equilibration.nf
 * Executes multi-stage energy minimization and equilibration (NVT/NPT) on Slurm GPU
 */

process GMX_EQUILIBRATION {
    tag "${meta.id}"
    label 'process_gpu'

    input:
    tuple val(meta), path(input_gro), path(topol_top), path(index_ndx), path(posre_files)
    path(mdp_minim)
    path(mdp_nvt)
    path(mdp_npt)
    val(is_membrane)
    val(temperature)

    output:
    tuple val(meta), path("eq_final.gro"), path("eq_final.cpt"), emit: equilibrated_state
    tuple val(meta), path("*.log"), emit: logs

    script:
    def pcoupl_type = is_membrane ? "semiisotropic" : "isotropic"
    def ref_p = is_membrane ? "1.0 1.0" : "1.0"
    def comp = is_membrane ? "4.5e-5 4.5e-5" : "4.5e-5"
    """
    echo "=== Worker: \$(hostname) - Starting Multi-stage Equilibration ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    # 1. Generate customized MDP files
    python3 - << 'EOF'
    import os
    is_mem = "${is_membrane}".lower() == "true"
    temp = "${temperature}"

    base_replace = {
        "{TEMP}": temp,
        "{TC_GRPS}": "PROT_LIG SOL_ION",
        "{PCOUPL_TYPE}": "${pcoupl_type}",
        "{REF_P}": "${ref_p}",
        "{COMPRESSIBILITY}": "${comp}"
    }

    if is_mem:
        stages = [
            ("emA", "${mdp_minim}", {"{DEFINE}": "-DPOSRES -DPOSRES_LIG", "{EMSTEP}": "0.005", "{CONSTRAINTS}": "h-bonds"}),
            ("emB", "${mdp_minim}", {"{DEFINE}": "", "{EMSTEP}": "0.001", "{CONSTRAINTS}": "h-bonds"}),
            ("emC", "${mdp_minim}", {"{DEFINE}": "", "{EMSTEP}": "0.0001", "{CONSTRAINTS}": "none"}),
            ("nvt_short", "${mdp_nvt}", {"{DEFINE}": "-DPOSRES -DPOSRES_LIG", "{NSTEPS}": "50000", "{DT}": "0.001", "{GEN_VEL}": "yes", "{CONTINUATION}": "no"}),
            ("nvt_long", "${mdp_nvt}", {"{DEFINE}": "-DPOSRES -DPOSRES_LIG", "{NSTEPS}": "250000", "{DT}": "0.002", "{GEN_VEL}": "no", "{CONTINUATION}": "yes"}),
            ("npt_berendsen", "${mdp_npt}", {"{DEFINE}": "-DPOSRES -DPOSRES_LIG", "{NSTEPS}": "500000", "{PCOUPL}": "C-rescale", "{TAU_P}": "5.0"}),
            ("npt_pr_250", "${mdp_npt}", {"{DEFINE}": "-DPOSRES_250 -DPOSRES_LIG_250", "{NSTEPS}": "500000", "{PCOUPL}": "Parrinello-Rahman", "{TAU_P}": "5.0"}),
            ("npt_pr_50", "${mdp_npt}", {"{DEFINE}": "-DPOSRES_50 -DPOSRES_LIG_50", "{NSTEPS}": "500000", "{PCOUPL}": "Parrinello-Rahman", "{TAU_P}": "5.0"}),
            ("npt_unrestrained", "${mdp_npt}", {"{DEFINE}": "", "{NSTEPS}": "500000", "{PCOUPL}": "Parrinello-Rahman", "{TAU_P}": "5.0"})
        ]
    else:
        stages = [
            ("emA", "${mdp_minim}", {"{DEFINE}": "-DPOSRES -DPOSRES_LIG", "{EMSTEP}": "0.01", "{CONSTRAINTS}": "h-bonds"}),
            ("nvt_short", "${mdp_nvt}", {"{DEFINE}": "-DPOSRES -DPOSRES_LIG", "{NSTEPS}": "50000", "{DT}": "0.001", "{GEN_VEL}": "yes", "{CONTINUATION}": "no"}),
            ("npt_pr_1000", "${mdp_npt}", {"{DEFINE}": "-DPOSRES -DPOSRES_LIG", "{NSTEPS}": "500000", "{PCOUPL}": "Parrinello-Rahman", "{TAU_P}": "5.0"}),
            ("npt_unrestrained", "${mdp_npt}", {"{DEFINE}": "", "{NSTEPS}": "500000", "{PCOUPL}": "Parrinello-Rahman", "{TAU_P}": "5.0"})
        ]

    for name, tpl, macros in stages:
        with open(tpl) as f: c = f.read()
        for k, v in base_replace.items(): c = c.replace(k, v)
        for k, v in macros.items(): c = c.replace(k, v)
        with open(f"{name}.mdp", "w") as out: out.write(c)

    with open("stages_list.txt", "w") as f:
        for name, _, _ in stages: f.write(f"{name}\n")
    EOF

    # 2. Sequential stage execution
    prev_gro="${input_gro}"
    prev_cpt=""

    while read stage; do
        echo "--> Running stage: \$stage"
        
        cpt_arg=""
        if [ -n "\$prev_cpt" ] && [ -f "\$prev_cpt" ] && [[ "\$stage" != *"em"* ]]; then
            cpt_arg="-t \$prev_cpt"
        fi

        gmx grompp -f "\${stage}.mdp" -o "\${stage}.tpr" -c "\$prev_gro" -r "\$prev_gro" -p "${topol_top}" -n "${index_ndx}" \$cpt_arg -maxwarn 2
        
        # Hardware offloading flags
        if [[ "\$stage" == *"em"* ]]; then
            gmx mdrun -deffnm "\$stage" -nb gpu -ntomp 8
        elif [[ "\$stage" == *"berendsen"* ]]; then
            gmx mdrun -deffnm "\$stage" -nb gpu -pme gpu -ntomp 8
        else
            gmx mdrun -deffnm "\$stage" -nb gpu -pme gpu -update gpu -ntomp 8
        fi

        prev_gro="\${stage}.gro"
        prev_cpt="\${stage}.cpt"
    done < stages_list.txt

    cp "\$prev_gro" eq_final.gro
    cp "\$prev_cpt" eq_final.cpt
    """

    stub:
    """
    touch eq_final.gro
    touch eq_final.cpt
    touch eq.log
    """
}

