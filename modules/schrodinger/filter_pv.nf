/*
 * modules/schrodinger/filter_pv.nf
 * Filters Pose Viewer file based on funnel whitelist CSV for a specific frame,
 * splitting into chunked PV files for massive distributed Prime MM-GBSA execution.
 * Self-contained worker process: zero external script file staging.
 */

process FILTER_PV {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(pv_maegz)
    path(whitelist_csv)
    val(batch_size)

    output:
    tuple val(meta), path("${meta.id}_filtered_pv.maegz"), emit: filtered_pv
    tuple val(meta), path("${meta.id}_filtered_c*.maegz"), emit: filtered_pv_chunks

    script:
    def frame_id = meta.frame ?: 1
    """
    cp \$(command -v filter_pv.py) .
    schrodinger run filter_pv.py \\
        --input ${pv_maegz} \\
        --output ${meta.id}_filtered_pv.maegz \\
        --out_prefix ${meta.id}_filtered \\
        --whitelist ${whitelist_csv} \\
        --frame ${frame_id} \\
        --batch_size ${batch_size} \\
        --max_poses 1
    """

    stub:
    """
    touch ${meta.id}_filtered_pv.maegz
    touch ${meta.id}_filtered_c0.maegz
    """
}
