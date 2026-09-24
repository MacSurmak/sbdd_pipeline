/*
 * modules/gromacs/run_chunk.nf
 * Orchestrates dynamic atomic production MD simulation chunks (1 to N chunks)
 * Completely eliminates any hardcoded chunk limits or manual DAG unrolling.
 */

process RUN_MD_CHUNK {
    tag "${meta.id}_R${meta.rep}_C${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), val(chunk_idx), path(prev_gro), path(prev_cpt), path(prev_tpr), path(topol_top), path(index_ndx)
    path mdp_templates
    val(is_membrane)
    val(temperature)
    val(chunk_ns)

    output:
    tuple val(meta), val(chunk_idx), path("md_${chunk_idx}.gro"), path("md_${chunk_idx}.cpt"), path("md_${chunk_idx}.tpr"), path(topol_top), path(index_ndx), emit: state
    tuple val(meta), path("md_${chunk_idx}.xtc"), emit: xtc
    tuple val(meta), path("sys_ref.tpr"), path(topol_top), path(index_ndx), emit: sys_files, optional: true

    script:
    def pcoupl_type = is_membrane.toString() == 'true' ? "semiisotropic" : "isotropic"
    def ref_p = is_membrane.toString() == 'true' ? "1.0 1.0" : "1.0"
    def comp = is_membrane.toString() == 'true' ? "4.5e-5 4.5e-5" : "4.5e-5"
    """
    echo "=== Worker: \$(hostname) | Chunk: ${chunk_idx} ==="
    
    CHUNK_STEPS=\$(awk "BEGIN {printf \\"%.0f\\", (${chunk_ns} * 1000.0) / 0.002}")
    EXTENSION_PS=\$(awk "BEGIN {printf \\"%.0f\\", ${chunk_ns} * 1000.0}")

    if [ "${chunk_idx}" -eq 1 ]; then
        cat md.mdp | sed -e "s/{NSTEPS}/\${CHUNK_STEPS}/g" \\
                         -e "s/{TEMP}/${temperature}/g" \\
                         -e "s/{TC_GRPS}/PROT_LIG SOL_ION/g" \\
                         -e "s/{PCOUPL_TYPE}/${pcoupl_type}/g" \\
                         -e "s/{REF_P}/${ref_p}/g" \\
                         -e "s/{COMPRESSIBILITY}/${comp}/g" > md_run.mdp

        gmx grompp -f md_run.mdp -o md_1.tpr -c ${prev_gro} -t ${prev_cpt} -p ${topol_top} -n ${index_ndx} -maxwarn 2
        gmx mdrun -deffnm md_1 -nb gpu -pme gpu -update gpu -ntomp ${task.cpus}
        cp md_1.tpr sys_ref.tpr
    else
        gmx convert-tpr -s ${prev_tpr} -extend \${EXTENSION_PS} -o md_${chunk_idx}.tpr
        gmx mdrun -deffnm md_${chunk_idx} -cpi ${prev_cpt} -noappend -nb gpu -pme gpu -update gpu -ntomp ${task.cpus}
    fi
    """

    stub:
    """
    touch md_${chunk_idx}.gro
    touch md_${chunk_idx}.cpt
    touch md_${chunk_idx}.tpr
    touch md_${chunk_idx}.xtc
    if [ "${chunk_idx}" -eq 1 ]; then
        touch sys_ref.tpr
    fi
    """
}
