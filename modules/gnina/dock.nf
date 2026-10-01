/*
 * modules/gnina/dock.nf
 * Distributed GNINA Docking Process for Slurm GPU nodes
 */

process GNINA_DOCK {
    tag "${meta.id}_chunk_${chunk_idx}"
    label 'process_gpu'

    input:
    tuple val(meta), path(receptor_pdb, stageAs: 'rec_dir/*'), path(ref_lig_sdf, stageAs: 'ref_dir/*'), path(chunk_sdf, stageAs: 'chunk_dir/*'), val(chunk_idx)
    val(exhaustiveness)
    val(cnn_scoring)
    val(num_modes)
    val(autobox_add)

    output:
    tuple val(meta), path("docked_chunk_${chunk_idx}.sdf"), emit: docked_sdf
    tuple val(meta), path("gnina_chunk_${chunk_idx}.log"), emit: log

    script:
    """
    REF_BASE=\$(basename "${ref_lig_sdf}")
    if [[ "\$REF_BASE" != "NO_FILE" && -f "${ref_lig_sdf}" ]]; then
        AUTOBOX_FLAG="--autobox_ligand ${ref_lig_sdf}"
    else
        AUTOBOX_FLAG=""
    fi

    echo "=== Worker: \$(hostname) (GPU Node) ==="
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

    gnina \
        -r ${receptor_pdb} \
        -l ${chunk_sdf} \
        \$AUTOBOX_FLAG \
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

