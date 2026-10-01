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
    tuple val(meta), path("${meta.id}_c*.pdb"), emit: cluster_pdbs

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

    # 2. Build index from dry reference PDB (Protein + Ligand + Cofactors)
    echo "q" | gmx make_ndx -f ${ref_pdb} -o index.ndx
    python3 ${projectDir}/bin/gromacs_ndx_fixer.py --ndx index.ndx --action ensure_complex

    # Fit on Protein scaffold, calculate RMSD matrix on selected group (default COMPLEX)
    printf "Protein\n${cluster_group}\n" | gmx cluster \\
        -f \${TARGET_XTC} \\
        -s ${ref_pdb} \\
        -n index.ndx \\
        -method ${method} \\
        -cutoff ${cutoff} \\
        ${dt_flag} \\
        -cl clusters_multi.pdb

    # 3. Split multi-model PDB into prefixed centroid PDBs
    python3 -c "
import sys
with open('clusters_multi.pdb', 'r') as f:
    lines = f.readlines()

out_f = None
idx = 1
for line in lines:
    if line.startswith('MODEL'):
        if out_f: out_f.close()
        out_f = open(f'${meta.id}_c{idx}.pdb', 'w')
        idx += 1
    if out_f:
        out_f.write(line)
if out_f: out_f.close()
"
    """

    stub:
    """
    touch ${meta.id}_c1.pdb
    touch ${meta.id}_c2.pdb
    """
}
