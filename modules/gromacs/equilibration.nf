/*
 * modules/gromacs/equilibration.nf
 * Executes multi-stage energy minimization and equilibration (NVT/NPT) on Slurm GPU
 */

process GMX_EQUILIBRATION {
    tag "${meta.id}"
    label 'process_gpu'

    input:
    tuple val(meta), path(input_gro), path(topol_top), path(index_ndx), path(posre_files)
    path md_protocol

    output:
    tuple val(meta), path("eq_final.gro"), path("eq_final.cpt"), path("eq_final.tpr"), emit: equilibrated_state
    tuple val(meta), path("*.log"), emit: logs

    script:
    """
    echo "=== Worker: \$(hostname) - Starting Multi-stage Equilibration ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    if [ ! -f protocol.txt ]; then
        echo "FATAL: protocol.txt playlist not found!"
        exit 1
    fi

    prev_gro="${input_gro}"
    prev_cpt=""
    prev_tpr=""

    while read mdp_file; do
        [ -z "\$mdp_file" ] && continue
        echo "--> Running stage: \$mdp_file"
        
        stage_name=\$(basename "\$mdp_file" .mdp)
        cpt_arg=""
        if [ -n "\$prev_cpt" ] && [ -f "\$prev_cpt" ] && [[ "\$stage_name" != *"em"* && "\$stage_name" != *"minim"* ]]; then
            cpt_arg="-t \$prev_cpt"
        fi

        gmx grompp -f "\$mdp_file" -o "\${stage_name}.tpr" -c "\$prev_gro" -r "\$prev_gro" -p "${topol_top}" -n "${index_ndx}" \$cpt_arg -maxwarn 2
        
        if [[ "\$stage_name" == *"em"* || "\$stage_name" == *"minim"* ]]; then
            gmx mdrun -deffnm "\${stage_name}" -nb gpu -ntmpi 1 -ntomp ${task.cpus}
        else
            gmx mdrun -deffnm "\${stage_name}" -nb gpu -pme gpu -ntmpi 1 -ntomp ${task.cpus}
        fi

        prev_gro="\${stage_name}.gro"
        prev_cpt="\${stage_name}.cpt"
        prev_tpr="\${stage_name}.tpr"
    done < protocol.txt

    cp "\$prev_gro" eq_final.gro
    cp "\$prev_cpt" eq_final.cpt
    cp "\$prev_tpr" eq_final.tpr
    """

    stub:
    """
    touch eq_final.gro
    touch eq_final.cpt
    touch eq_final.tpr
    touch eq.log
    """
}
