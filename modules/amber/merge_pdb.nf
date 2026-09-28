/*
 * modules/amber/merge_pdb.nf
 * Cleanly merges protein and ligand PDB files using native bash utilities.
 */

process MERGE_PDB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_pdb), path(org_pdbs)

    output:
    tuple val(meta), path("complex_dry.pdb"), emit: complex_pdb

    script:
    """
    grep -v -E "^(END|TER)" ${rec_pdb} > complex_dry.pdb
    for f in ${org_pdbs}; do
        grep -E "^(ATOM|HETATM)" \$f >> complex_dry.pdb
    done
    echo "TER" >> complex_dry.pdb
    echo "END" >> complex_dry.pdb
    """

    stub:
    """
    touch complex_dry.pdb
    """
}
