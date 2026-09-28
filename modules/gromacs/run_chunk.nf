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
    path md_protocol
    val(chunk_ns)
    val(chunk_idx)

    output:
    tuple val(meta), path("out.gro"), path("out.cpt"), path("out.tpr"), path(topol_top), path(index_ndx), emit: state
    tuple val(meta), path("md_${chunk_idx}.xtc"), emit: xtc
    tuple val(meta), path("sys_ref.tpr"), path(topol_top), path(index_ndx), emit: sys_files

    script:
    """
    echo "=== Worker: \$(hostname) | Chunk: ${chunk_idx} ==="
    
    CHUNK_STEPS=\$(awk "BEGIN {printf \\"%.0f\\", (${chunk_ns} * 1000.0) / 0.002}")
    EXTENSION_PS=\$(awk "BEGIN {printf \\"%.0f\\", ${chunk_ns} * 1000.0}")
    NSTXOUT=\$(awk "BEGIN {steps=\${CHUNK_STEPS}; n=int(steps/10); if (n<1) n=1; if (n>5000) n=5000; print n}")
    NSTENERGY=\$(awk "BEGIN {steps=\${CHUNK_STEPS}; n=int(steps/10); if (n<1) n=1; if (n>1000) n=1000; print n}")

    if [ "${chunk_idx}" -eq 1 ]; then
        sed -e "s/nsteps\\s*=.*/nsteps = \${CHUNK_STEPS}/g" \\
            -e "s/nstxout-compressed\\s*=.*/nstxout-compressed = \${NSTXOUT}/g" \\
            -e "s/nstenergy\\s*=.*/nstenergy = \${NSTENERGY}/g" \\
            md.mdp > md_run.mdp

        gmx grompp -f md_run.mdp -o out.tpr -c ${prev_gro} -t ${prev_cpt} -p ${topol_top} -n ${index_ndx} -maxwarn 2
        gmx mdrun -deffnm md_${chunk_idx} -s out.tpr -nb gpu -pme gpu -ntmpi 1 -ntomp ${task.cpus}
        cp out.tpr sys_ref.tpr
    else
        gmx convert-tpr -s ${prev_tpr} -extend \${EXTENSION_PS} -o out.tpr
        gmx mdrun -deffnm md_${chunk_idx} -s out.tpr -cpi ${prev_cpt} -noappend -nb gpu -pme gpu -ntmpi 1 -ntomp ${task.cpus}
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
