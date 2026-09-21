/*
 * modules/gromacs/run_chunk.nf
 * Executes a single production MD chunk (e.g. 10 ns) on Slurm GPU
 */

process RUN_MD_CHUNK {
    tag "${meta.id}_rep${meta.rep}_chunk${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), path(in_gro), path(in_cpt), path(in_tpr), path(topol_top), path(index_ndx), val(chunk_idx)
    path(mdp_template)
    val(chunk_ns)
    val(temperature)

    output:
    tuple val(meta), path("md_${chunk_idx}.gro"), path("md_${chunk_idx}.cpt"), path("md_${chunk_idx}.tpr"), val(chunk_idx), emit: next_chunk_input
    tuple val(meta), path("md_${chunk_idx}.xtc"), emit: xtc
    tuple val(meta), path("md_${chunk_idx}.edr"), emit: edr
    tuple val(meta), path("md_${chunk_idx}.log"), emit: log

    script:
    def extend_ps = (chunk_ns * 1000) as int
    def nsteps = ((chunk_ns * 1000) / 0.002) as int
    """
    echo "=== Worker: \$(hostname) - Executing Chunk ${chunk_idx} for ${meta.id} (Rep ${meta.rep}) ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    if [ "${chunk_idx}" -eq 1 ]; then
        # Prepare customized md.mdp with proper steps and temperature
        sed -e 's/{TEMP}/${temperature}/g' \
            -e 's/{TC_GRPS}/PROT_LIG SOL_ION/g' \
            -e 's/{NSTEPS}/${nsteps}/g' \
            ${mdp_template} > md.mdp

        # Seed randomization for independent replicates
        echo "gen-vel = yes" >> md.mdp
        echo "gen-temp = ${temperature}" >> md.mdp
        echo "gen-seed = ${meta.rep * 12345 + 101}" >> md.mdp

        gmx grompp -f md.mdp -o md_1.tpr -c ${in_gro} -t ${in_cpt} -p ${topol_top} -n ${index_ndx} -maxwarn 2
        gmx mdrun -deffnm md_1 -nb gpu -pme gpu -update gpu -ntomp 8 -noappend
    else
        gmx convert-tpr -s ${in_tpr} -extend ${extend_ps} -o md_${chunk_idx}.tpr
        gmx mdrun -deffnm md_${chunk_idx} -cpi ${in_cpt} -nb gpu -pme gpu -update gpu -ntomp 8 -noappend
    fi
    """

    stub:
    """
    touch md_${chunk_idx}.gro
    touch md_${chunk_idx}.cpt
    touch md_${chunk_idx}.tpr
    touch md_${chunk_idx}.xtc
    touch md_${chunk_idx}.edr
    touch md_${chunk_idx}.log
    """
}

