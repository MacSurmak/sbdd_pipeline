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

    schrodinger run python3 -c "
import os
import sys
from pathlib import Path
from schrodinger import structure

lig_name = '${meta.ligand}'.strip().upper()
max_c = int('${meta.max_clusters}' or 0)
receipt_lines = []

with structure.StructureReader('clusters_multi.pdb') as reader:
    for i, st in enumerate(reader, 1):
        if max_c > 0 and i > max_c:
            break
            
        pose_path = 'NO_POSE'
        if lig_name:
            lig_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == lig_name]
            if lig_atoms:
                lig_st = st.extract(lig_atoms)
                lig_st.title = lig_name
                lig_st.property['s_m_title'] = lig_name
                sdf_out = f'${meta.id}_pose_{i}.sdf'
                mae_out = f'${meta.id}_pose_{i}.mae'
                lig_st.write(sdf_out)
                lig_st.write(mae_out)
                st.deleteAtoms(lig_atoms)
                pose_path = str(Path(sdf_out).resolve())
                
        c_title = f'${meta.id}_centroid_{i}'
        st.title = c_title
        st.property['s_m_title'] = c_title
        pdb_out = f'${meta.id}_centroid_{i}.pdb'
        mae_out = f'${meta.id}_centroid_{i}.mae'
        st.write(pdb_out)
        st.write(mae_out)
        
        receipt_lines.append(f'{c_title}\t{Path(pdb_out).resolve()}\t{Path(mae_out).resolve()}\t{pose_path}\n')

with open('cluster_receipt_chunk.tsv', 'w') as rf:
    rf.writelines(receipt_lines)
"
    """

    stub:
    """
    touch ${meta.id}_centroid_1.pdb
    touch ${meta.id}_centroid_1.mae
    touch ${meta.id}_pose_1.sdf
    touch ${meta.id}_pose_1.mae
    echo -e "${meta.id}_centroid_1\t$PWD/${meta.id}_centroid_1.pdb\t$PWD/${meta.id}_centroid_1.mae\t$PWD/${meta.id}_pose_1.sdf" > cluster_receipt_chunk.tsv
    """
}
