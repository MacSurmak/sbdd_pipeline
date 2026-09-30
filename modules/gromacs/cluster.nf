/*
 * modules/gromacs/cluster.nf
 * GROMACS trajectory clustering and automatic PDB splitting
 */

process GMX_CLUSTER {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(xtc), path(tpr), path(ref_pdb)
    val(cutoff)
    val(method)

    output:
    tuple val(meta), path("cluster_*.pdb"), emit: cluster_pdbs

    script:
    """
    # 1. Create strict index with Protein group if needed
    echo "q" | gmx make_ndx -f ${ref_pdb} -o index.ndx

    # 2. Run clustering (fit on Protein, calculate RMSD on Protein)
    echo "Protein Protein" | gmx cluster -f ${xtc} -s ${tpr} -n index.ndx -method ${method} -cutoff ${cutoff} -cl clusters_multi.pdb

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
