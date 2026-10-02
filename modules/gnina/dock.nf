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
    tuple val(meta), path("${meta.id}_docked.sdf"), emit: docked_sdf
    tuple val(meta), path("${meta.id}_gnina.log"), emit: log

    script:
    def autobox_flag = ""
    """
    if [ ! -f "${ref_lig_sdf}" ] || [ ! -s "${ref_lig_sdf}" ]; then
        echo "FATAL: Reference ligand SDF is missing or empty: ${ref_lig_sdf}"
        exit 1
    fi
    AUTOBOX_FLAG="--autobox_ligand ${ref_lig_sdf}"

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
        -o ${meta.id}_docked.sdf > ${meta.id}_gnina.log
    """

    stub:
    """
    touch ${meta.id}_docked.sdf
    touch ${meta.id}_gnina.log
    """}

