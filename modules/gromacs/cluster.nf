/*
 * modules/gromacs/cluster.nf
 * GROMACS trajectory clustering and automatic PDB splitting
 */

process GMX_CLUSTER {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(clean_xtcs, stageAs: 'trajs_?/*'), path(ref_pdb)
    val(cutoff)
    val(method)

    output:
    tuple val(meta), path("cluster_*.pdb"), emit: cluster_pdbs

    script:
    """
    # 1. Discover clean, pre-centered, solvent-stripped trajectories
    XTC_FILES=(trajs_*/*.xtc)
    if [ \${#XTC_FILES[@]} -gt 1 ]; then
        echo "Pooling \${#XTC_FILES[@]} clean replica trajectories into single conformational ensemble..."
        gmx trjcat -f \${XTC_FILES[@]} -o combined_clean.xtc -cat
        TARGET_XTC="combined_clean.xtc"
    else
        TARGET_XTC="\${XTC_FILES[0]}"
    fi

    # 2. Cluster the complete COMPLEX (Protein + Ligand + Cofactors)
    # Generates strict index from the dry reference PDB
    echo "q" | gmx make_ndx -f ${ref_pdb} -o index.ndx

    # Fit on Protein, calculate RMSD on Protein, output the full COMPLEX into cluster centroids
    printf "Protein\nProtein\n" | gmx cluster -f \${TARGET_XTC} -s ${ref_pdb} -n index.ndx -method ${method} -cutoff ${cutoff} -cl clusters_multi.pdb

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
