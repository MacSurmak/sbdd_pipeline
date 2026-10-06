process GENERATE_DOCKING_BOX {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(ref_sdf)
    val(inner_box)
    val(outer_box)
    val(autobox_add)

    output:
    tuple val(meta), path("${meta.id}_docking_box.pml"), emit: box_pml

    script:
    """
    python3 ${projectDir}/bin/generate_docking_box.py \\
        --reference ${ref_sdf} \\
        --receptor_name "${meta.id}" \\
        --inner_box "${inner_box}" \\
        --outer_box "${outer_box}" \\
        --autobox_add ${autobox_add} \\
        --output "${meta.id}_docking_box.pml"
    """
}

process CALC_RMSD {
    tag "${meta.id}_${engine}"
    label 'process_single'

    input:
    tuple val(meta), path(docked_sdf), path(ref_sdf)
    val(engine)

    output:
    tuple val(meta), path("${meta.id}_${engine}_rmsd.csv"), emit: rmsd_csv

    script:
    def lig_title = meta.native_ligand ?: meta.ligand_title ?: ""
    """
    python3 ${projectDir}/bin/calc_rmsd.py \\
        --reference ${ref_sdf} \\
        --poses ${docked_sdf} \\
        --receptor "${meta.frame_id ?: meta.id}" \\
        --title "${lig_title}" \\
        --engine ${engine} \\
        --output ${meta.id}_${engine}_rmsd.csv
    """
}

process MERGE_GNINA_DOCK {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(chunk_sdfs)

    output:
    tuple val(meta), path("${meta.id}_gnina_docked_all.sdf"), emit: merged_sdf

    script:
    """
    python3 ${projectDir}/bin/wrap_utils.py merge_sdf --in_sdfs ${chunk_sdfs} --out_sdf ${meta.id}_gnina_docked_all.sdf
    """
}

process MERGE_GNINA_RESC {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(chunk_sdfs)

    output:
    tuple val(meta), path("${meta.id}_rescored.sdf"), emit: rescored_sdf

    script:
    """
    python3 ${projectDir}/bin/wrap_utils.py merge_sdf --in_sdfs ${chunk_sdfs} --out_sdf ${meta.id}_rescored.sdf
    """
}

process MERGE_GLIDE {
    tag "${meta.id}"
    executor 'local'
    scratch false
    label 'process_single'

    input:
    tuple val(meta), path(pv_files, stageAs: 'pvs_?/*'), path(csv_files, stageAs: 'csvs_?/*'), path(sdf_files, stageAs: 'sdfs_?/*')

    output:
    tuple val(meta), path("${meta.id}_dock_pv.maegz"), emit: pv_maegz
    tuple val(meta), path("${meta.id}_dock.csv"), emit: scores_csv
    tuple val(meta), path("${meta.id}_dock.sdf"), emit: poses_sdf

    script:
    """
    schrodinger run ${projectDir}/bin/merge_pv.py \\
        --pv_files ${pv_files} \\
        --csv_files ${csv_files} \\
        --sdf_files ${sdf_files} \\
        --out_pv ${meta.id}_dock_pv.maegz \\
        --out_csv ${meta.id}_dock.csv \\
        --out_sdf ${meta.id}_dock.sdf
    """
}
