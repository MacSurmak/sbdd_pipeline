/*
 * modules/amber/merge_pdb.nf
 * Cleanly merges protein and ligand PDB files using native bash utilities.
 */

process MERGE_PDB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_pdb), path(lig_sdf)

    output:
    tuple val(meta), path("complex_dry.pdb"), emit: complex_pdb

    script:
    def lig_arg = lig_sdf ? "--lig_sdf ${lig_sdf}" : ""
    """
    python3 ${projectDir}/bin/merge_complex.py \\
        --rec_pdb ${rec_pdb} \\
        --out_pdb complex_dry.pdb \\
        ${lig_arg}
    """

    stub:
    """
    touch complex_dry.pdb
    """
}
