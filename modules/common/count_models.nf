/*
 * modules/common/count_models.nf
 * Automatically counts structures or NMR models in macromolecular file (.cif / .pdb / .mae)
 */

process COUNT_MODELS {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(receptor_file)

    output:
    tuple val(meta), path("models_count.txt"), emit: count_file

    script:
    """
    ${projectDir}/bin/count_models.py \\
        --input ${receptor_file} \\
        --output models_count.txt
    """

    stub:
    """
    echo "2" > models_count.txt
    """
}
