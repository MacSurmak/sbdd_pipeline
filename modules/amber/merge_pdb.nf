/*
 * modules/amber/merge_pdb.nf
 * Cleanly merges protein and ligand PDB files using native bash utilities.
 */

process MERGE_PDB {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_pdb), path(lig_sdf), path(extra_sdfs)

    output:
    tuple val(meta), path("complex_dry.pdb"), emit: complex_pdb

    script:
    def lig_arg = (lig_sdf && lig_sdf.name != 'NO_FILE') ? "--lig_sdf ${lig_sdf} --lig_resname ${meta.lig_resname}" : ""
    def extra_sdfs_arg = (extra_sdfs && extra_sdfs.size() > 0) ? "--extra_sdfs ${extra_sdfs}" : ""
    def extra_res_arg = (meta.extra_resnames && meta.extra_resnames.size() > 0) ? "--extra_resnames ${meta.extra_resnames.join(' ')}" : ""
    """
    python3 ${projectDir}/bin/merge_complex.py \\
        --rec_pdb ${rec_pdb} \\
        --out_pdb complex_dry.pdb \\
        ${lig_arg} \\
        ${extra_sdfs_arg} \\
        ${extra_res_arg}
    """
}
