/*
 * modules/gromacs/cluster.nf
 * GROMACS trajectory clustering and automatic PDB splitting
 */

process GMX_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(xtcs), path(tpr)
    val(cutoff)
    val(method)

    output:
    tuple val(meta), path("cluster_*.pdb"), emit: cluster_pdbs

    script:
    """
    # 1. Generate index strictly from binary system topology (TPR)
    echo "q" | gmx make_ndx -f ${tpr} -o index.ndx

    # 2. Automatically concatenate multiple replicas if provided
    XTC_ARRAY=(${xtcs})
    if [ \${#XTC_ARRAY[@]} -gt 1 ]; then
        echo "Concatenating \${#XTC_ARRAY[@]} replica trajectories for ensemble clustering..."
        gmx trjcat -f ${xtcs} -o combined_ensemble.xtc -cat
        TARGET_XTC="combined_ensemble.xtc"
    else
        TARGET_XTC="${xtcs}"
    fi

    # 3. Run clustering (fit on Protein, RMSD on Protein, output Protein to clusters)
    printf "Protein\nProtein\n" | gmx cluster -f \${TARGET_XTC} -s ${tpr} -n index.ndx -method ${method} -cutoff ${cutoff} -cl clusters_multi.pdb

    # 3. Split multi-model PDB into separate files for database registration
    python3 -c "
import sys
with open('clusters_multi.pdb', 'r') as f:
    lines = f.readlines()

out_f = None
idx = 1
for line in lines:
    if line.startswith('MODEL'):
        if out_f: out_f.close()
        out_f = open(f'cluster_{idx}.pdb', 'w')
        idx += 1
    if out_f:
        out_f.write(line)
if out_f: out_f.close()
"
    """

    stub:
    """
    touch cluster_1.pdb
    touch cluster_2.pdb
    """
}
