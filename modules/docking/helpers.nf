process ASSEMBLE_DOCKING_RECEPTOR {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_pdb), path(rec_mae), path(ion_sdfs), path(ion_maes), path(cof_sdfs), path(cof_maes)

    output:
    tuple val(meta), path("${meta.id}_dock_receptor.pdb"), path("${meta.id}_dock_receptor.mae"), emit: assembled_receptors

    script:
    def ion_sdfs_arg = (ion_sdfs && ion_sdfs.size() > 0) ? "--ion_sdfs ${ion_sdfs}" : ""
    def ion_res_arg = (meta.ion_resnames && meta.ion_resnames.size() > 0) ? "--ion_resnames ${meta.ion_resnames.join(' ')}" : ""
    def cof_sdfs_arg = (cof_sdfs && cof_sdfs.size() > 0) ? "--cofactor_sdfs ${cof_sdfs}" : ""
    def cof_res_arg = (meta.cofactor_resnames && meta.cofactor_resnames.size() > 0) ? "--cofactor_resnames ${meta.cofactor_resnames.join(' ')}" : ""

    def all_hetero_maes = []
    if (ion_maes && ion_maes.size() > 0) all_hetero_maes.addAll(ion_maes instanceof List ? ion_maes : [ion_maes])
    if (cof_maes && cof_maes.size() > 0) all_hetero_maes.addAll(cof_maes instanceof List ? cof_maes : [cof_maes])
    def hetero_maes_arg = all_hetero_maes ? "--hetero_maes ${all_hetero_maes.join(' ')}" : ""

    """
    # PDB Assembly (for GNINA and scoring)
    python3 ${projectDir}/bin/merge_complex.py \\
        --rec_pdb ${rec_pdb} \\
        --out_pdb ${meta.id}_dock_receptor.pdb \\
        ${ion_sdfs_arg} \\
        ${ion_res_arg} \\
        ${cof_sdfs_arg} \\
        ${cof_res_arg}

    # MAE Assembly (for Glide and PrepWizard state preservation)
    schrodinger run python3 ${projectDir}/bin/assemble_mae.py \\
        --rec_mae ${rec_mae} \\
        --out_mae ${meta.id}_dock_receptor.mae \\
        ${hetero_maes_arg}
    """
}

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
