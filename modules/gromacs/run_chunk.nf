/*
 * modules/gromacs/run_chunk.nf
 * Orchestrates dynamic atomic production MD simulation chunks (1 to N chunks)
 * Completely eliminates any hardcoded chunk limits or manual DAG unrolling.
 */

process RUN_MD_CHUNK {
    tag "${meta.id}_R${meta.rep}_C${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), path(prev_gro), path(prev_cpt), path(prev_tpr), path(topol_top), path(index_ndx)
    path mdp_templates
    val(is_membrane)
    val(temperature)
    val(chunk_ns)
    val(chunk_idx)

    output:
    tuple val(meta), path("out.gro"), path("out.cpt"), path("out.tpr"), path(topol_top), path(index_ndx), emit: state
    tuple val(meta), path("md_${chunk_idx}.xtc"), emit: xtc
    tuple val(meta), path("sys_ref.tpr"), path(topol_top), path(index_ndx), emit: sys_files

    script:
    def pcoupl_type = is_membrane.toString() == 'true' ? "semiisotropic" : "isotropic"
    def ref_p = is_membrane.toString() == 'true' ? "1.0 1.0" : "1.0"
    def comp = is_membrane.toString() == 'true' ? "4.5e-5 4.5e-5" : "4.5e-5"
    """
    echo "=== Worker: \$(hostname) | Chunk: ${chunk_idx} ==="
    
    CHUNK_STEPS=\$(awk "BEGIN {printf \\"%.0f\\", (${chunk_ns} * 1000.0) / 0.002}")
    EXTENSION_PS=\$(awk "BEGIN {printf \\"%.0f\\", ${chunk_ns} * 1000.0}")

    if [ "${chunk_idx}" -eq 1 ]; then
        cat md.mdp | sed -e "s/{NSTEPS}/\${CHUNK_STEPS}/g" \
                         -e "s/{TEMP}/${temperature}/g" \
                         -e "s/{TC_GRPS}/PROT_LIG SOL_ION/g" \
                         -e "s/{PCOUPL_TYPE}/${pcoupl_type}/g" \
                         -e "s/{REF_P}/${ref_p}/g" \
                         -e "s/{COMPRESSIBILITY}/${comp}/g" > md_run.mdp

        gmx grompp -f md_run.mdp -o out.tpr -c ${prev_gro} -t ${prev_cpt} -p ${topol_top} -n ${index_ndx} -maxwarn 2
        gmx mdrun -deffnm md_${chunk_idx} -s out.tpr -nb gpu -pme gpu -ntmpi 1 -ntomp ${task.cpus}
        cp out.tpr sys_ref.tpr
    else
        gmx convert-tpr -s ${prev_tpr} -extend \${EXTENSION_PS} -o out.tpr
        gmx mdrun -deffnm md_${chunk_idx} -s out.tpr -cpi ${prev_cpt} -nb gpu -pme gpu -ntmpi 1 -ntomp ${task.cpus}
        cp out.tpr sys_ref.tpr
    fi

    if [ -f md_${chunk_idx}.part*.gro ]; then
        mv md_${chunk_idx}.part*.gro md_${chunk_idx}.gro
    fi
    if [ -f md_${chunk_idx}.part*.cpt ]; then
        mv md_${chunk_idx}.part*.cpt md_${chunk_idx}.cpt
    fi
    if [ -f md_${chunk_idx}.part*.xtc ]; then
        mv md_${chunk_idx}.part*.xtc md_${chunk_idx}.xtc
    fi

    cp md_${chunk_idx}.gro out.gro
    cp md_${chunk_idx}.cpt out.cpt
    """

    stub:
    """
    touch out.gro
    touch out.cpt
    touch out.tpr
    touch md_${chunk_idx}.xtc
    touch sys_ref.tpr
    """
}
