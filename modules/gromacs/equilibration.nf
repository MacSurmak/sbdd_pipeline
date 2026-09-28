/*
 * modules/gromacs/equilibration.nf
 * Executes multi-stage energy minimization and equilibration (NVT/NPT) on Slurm GPU
 */

process GMX_EQUILIBRATION {
    tag "${meta.id}"
    label 'process_gpu'

    input:
    tuple val(meta), path(input_gro), path(topol_top), path(index_ndx), path(posre_files)
    path mdp_templates
    val(is_membrane)
    val(temperature)
    val(nvt_short_ps)
    val(nvt_long_ps)
    val(npt_berendsen_ps)
    val(npt_pr_ps)

    output:
    tuple val(meta), path("eq_final.gro"), path("eq_final.cpt"), emit: equilibrated_state
    tuple val(meta), path("*.log"), emit: logs

    script:
    """
    echo "=== Worker: \$(hostname) - Starting Multi-stage Equilibration ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    # 1. Generate customized MDP files from templates
    python3 ${projectDir}/bin/generate_mdp.py \\
        --templates_dir . \\
        --outdir . \\
        --is_membrane ${is_membrane} \\
        --temperature ${temperature} \\
        --nvt_short_ps ${nvt_short_ps} \\
        --nvt_long_ps ${nvt_long_ps} \\
        --npt_berendsen_ps ${npt_berendsen_ps} \\
        --npt_pr_ps ${npt_pr_ps}

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
        
        # Hardware offloading flags (1 thread-MPI rank per GPU + task.cpus OpenMP threads)
        if [[ "\$stage" == *"em"* ]]; then
            gmx mdrun -deffnm "\$stage" -nb gpu -ntmpi 1 -ntomp ${task.cpus}
        else
            gmx mdrun -deffnm "\$stage" -nb gpu -pme gpu -ntmpi 1 -ntomp ${task.cpus}
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
