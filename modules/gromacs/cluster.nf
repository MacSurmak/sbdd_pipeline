/*
 * modules/gromacs/cluster.nf
 * GROMACS trajectory clustering and automatic PDB splitting
 */

process GMX_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(clean_xtcs, stageAs: 'trajs_?/*'), path(ref_pdb), val(cutoff), val(method), val(dt), val(cluster_group)

    output:
    tuple val(meta), path("${meta.id}_centroid_*.pdb"), path("${meta.id}_centroid_*.mae"), emit: cluster_frames
    tuple val(meta), path("${meta.id}_pose_*.sdf"), path("${meta.id}_pose_*.mae"), emit: cluster_poses, optional: true
    path("cluster_receipt_chunk.tsv"), emit: cluster_manifest

    script:
    def dt_flag = (dt && dt.toString() != '0') ? "-dt ${dt}" : ""
    """
    XTC_FILES=(trajs_*/*.xtc)
    if [ \${#XTC_FILES[@]} -gt 1 ]; then
        echo "Pooling \${#XTC_FILES[@]} clean replica trajectories into single conformational ensemble..."
        gmx trjcat -f \${XTC_FILES[@]} -o combined_clean.xtc -cat
        TARGET_XTC="combined_clean.xtc"
    else
        TARGET_XTC="\${XTC_FILES[0]}"
    fi

    echo "q" | gmx make_ndx -f ${ref_pdb} -o index.ndx
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action ensure_complex

    TARGET_GROUP="${cluster_group}"
    if [ "${cluster_group}" = "Pocket" ] || [ "${cluster_group}" = "pocket" ]; then
        python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --pdb ${ref_pdb} --ligand "${meta.ligand}" --action add_pocket
        TARGET_GROUP="Pocket"
    elif [ "${cluster_group}" = "Ligand" ] || [ "${cluster_group}" = "ligand" ]; then
        TARGET_GROUP="${meta.ligand}"
    fi

    wcl_flag=""
    if [ -n "${meta.max_clusters}" ] && [ "${meta.max_clusters}" -gt 0 ]; then
        wcl_flag="-wcl ${meta.max_clusters}"
    fi

    printf "\${TARGET_GROUP}\nCOMPLEX\n" | gmx cluster \
        -f \${TARGET_XTC} \
        -s ${ref_pdb} \
        -n index.ndx \
        -method ${method} \
        -cutoff ${cutoff} \
        ${dt_flag} \
        \${wcl_flag} \
        -cl clusters_multi.pdb

    schrodinger run python3 ${projectDir}/bin/split_cluster_pdb.py \
        --input clusters_multi.pdb \
        --prefix "${meta.id}" \
        --ligand "${meta.ligand ?: ''}" \
        --max_clusters "${meta.max_clusters ?: 0}" \
        --out_receipt cluster_receipt_chunk.tsv
    """

    stub:
    """
    touch ${meta.id}_centroid_1.pdb
    touch ${meta.id}_centroid_1.mae
    touch ${meta.id}_pose_1.sdf
    touch ${meta.id}_pose_1.mae
    printf "%s\\t%s/%s\\t%s/%s\\t%s/%s\\n" "${meta.id}_centroid_1" "\$PWD" "${meta.id}_centroid_1.pdb" "\$PWD" "${meta.id}_centroid_1.mae" "\$PWD" "${meta.id}_pose_1.sdf" > cluster_receipt_chunk.tsv
    """
}
