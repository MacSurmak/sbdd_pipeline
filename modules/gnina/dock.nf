/*
 * modules/gnina/dock.nf
 * Distributed GNINA Docking Process for Slurm GPU nodes
 */

process GNINA_DOCK {
    tag "${meta.id}_chunk_${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), path(receptor_pdb), path(ref_lig_sdf), path(chunk_sdf), val(chunk_idx)
    val(exhaustiveness)
    val(cnn_scoring)
    val(num_modes)
    val(autobox_add)

    output:
    tuple val(meta), path("docked_chunk_${chunk_idx}.sdf"), emit: docked_sdf
    tuple val(meta), path("gnina_chunk_${chunk_idx}.log"), emit: log

    script:
    """
    echo "=== Worker: \$(hostname) (GPU Node) ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    gnina \
        -r ${receptor_pdb} \
        -l ${chunk_sdf} \
        --autobox_ligand ${ref_lig_sdf} \
        --autobox_add ${autobox_add} \
        --exhaustiveness ${exhaustiveness} \
        --cnn_scoring ${cnn_scoring} \
        --num_modes ${num_modes} \
        -o docked_chunk_${chunk_idx}.sdf > gnina_chunk_${chunk_idx}.log
    """

    stub:
    """
    touch docked_chunk_${chunk_idx}.sdf
    touch gnina_chunk_${chunk_idx}.log
    """
}

