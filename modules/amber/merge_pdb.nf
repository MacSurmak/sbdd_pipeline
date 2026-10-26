/*
 * modules/amber/merge_pdb.nf
 * Cleanly merges protein and ligand PDB files using native bash utilities.
 */

process MERGE_PDB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_pdb), path(lig_sdf), path(ion_sdfs), path(cofactor_sdfs)

    output:
    tuple val(meta), path("complex_dry.pdb"), emit: complex_pdb

    script:
    def lig_arg = (lig_sdf && lig_sdf.name != 'NO_FILE') ? "--lig_sdf ${lig_sdf} --lig_resname ${meta.lig_resname}" : ""
    def ion_sdfs_arg = (ion_sdfs && ion_sdfs.size() > 0) ? "--ion_sdfs ${ion_sdfs}" : ""
    def ion_res_arg = (meta.ion_resnames && meta.ion_resnames.size() > 0) ? "--ion_resnames ${meta.ion_resnames.join(' ')}" : ""
    def cof_sdfs_arg = (cofactor_sdfs && cofactor_sdfs.size() > 0) ? "--cofactor_sdfs ${cofactor_sdfs}" : ""
    def cof_res_arg = (meta.cofactor_resnames && meta.cofactor_resnames.size() > 0) ? "--cofactor_resnames ${meta.cofactor_resnames.join(' ')}" : ""
    """
    python3 ${projectDir}/bin/merge_complex.py \\
        --rec_pdb ${rec_pdb} \\
        --out_pdb complex_dry.pdb \\
        ${lig_arg} \\
        ${ion_sdfs_arg} \\
        ${ion_res_arg} \\
        ${cof_sdfs_arg} \\
        ${cof_res_arg}
    """
}
