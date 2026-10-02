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

    script:
    def dt_flag = (dt && dt.toString() != '0') ? "-dt ${dt}" : ""
    """
    # 1. Discover and pool clean replica trajectories into a single conformational ensemble
    XTC_FILES=(trajs_*/*.xtc)
    if [ \${#XTC_FILES[@]} -gt 1 ]; then
        echo "Pooling \${#XTC_FILES[@]} clean replica trajectories into single conformational ensemble..."
        gmx trjcat -f \${XTC_FILES[@]} -o combined_clean.xtc -cat
        TARGET_XTC="combined_clean.xtc"
    else
        TARGET_XTC="\${XTC_FILES[0]}"
    fi

    # 2. Build index and dynamic Pocket group if requested
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

    # Fit and calculate RMSD strictly on TARGET_GROUP, output full COMPLEX to PDB
    printf "\${TARGET_GROUP}\nCOMPLEX\n" | gmx cluster \
        -f \${TARGET_XTC} \
        -s ${ref_pdb} \
        -n index.ndx \
        -method ${method} \
        -cutoff ${cutoff} \
        ${dt_flag} \
        \${wcl_flag} \
        -cl clusters_multi.pdb

    # 3. Split multi-model PDB and generate proper MAE/SDF formats with Schrodinger API
    schrodinger run python3 -c "
import sys
from schrodinger import structure

lig_name = '${meta.ligand}'.strip().upper()
max_c = int('${meta.max_clusters}' or 0)

with structure.StructureReader('clusters_multi.pdb') as reader:
    for i, st in enumerate(reader, 1):
        if max_c > 0 and i > max_c:
            break
            
        if lig_name:
            lig_atoms = [a.index for a in st.atom if a.pdbres.strip().upper() == lig_name]
            if lig_atoms:
                lig_st = st.extract(lig_atoms)
                lig_st.title = lig_name
                lig_st.property['s_m_title'] = lig_name
                lig_st.write(f'${meta.id}_pose_{i}.sdf')
                lig_st.write(f'${meta.id}_pose_{i}.mae')
                st.deleteAtoms(lig_atoms)
                
        st.title = f'${meta.id}_centroid_{i}'
        st.property['s_m_title'] = f'${meta.id}_centroid_{i}'
        st.write(f'${meta.id}_centroid_{i}.pdb')
        st.write(f'${meta.id}_centroid_{i}.mae')
"
    """

    stub:
    """
    touch ${meta.id}_centroid_1.pdb
    touch ${meta.id}_centroid_1.mae
    touch ${meta.id}_pose_1.sdf
    touch ${meta.id}_pose_1.mae
    """}
