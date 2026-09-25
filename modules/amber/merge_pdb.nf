/*
 * modules/amber/merge_pdb.nf
 * Cleanly merges protein and ligand PDB files using native bash utilities.
 */

process MERGE_PDB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_amber_pdb), path(unk_pdb)

    output:
    tuple val(meta), path("complex_dry.pdb"), emit: complex_pdb

    script:
    """
    # Strip trailing END/TER records from the protein
    grep -v -E "^(END|TER)" ${rec_amber_pdb} > complex_dry.pdb
    
    # Append purely structural records from the ligand
    grep -E "^(ATOM|HETATM)" ${unk_pdb} >> complex_dry.pdb
    
    # Re-seal the PDB format
    echo "TER" >> complex_dry.pdb
    echo "END" >> complex_dry.pdb
    """

    stub:
    """
    touch complex_dry.pdb
    """
}
