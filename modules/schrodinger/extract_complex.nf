/*
 * modules/schrodinger/extract_complex.nf
 * Extracts raw protein PDB and native ligand MAE & SDF from complex (.cif / .pdb)
 * Self-contained worker process: zero external script file staging.
 */

process EXTRACT_COMPLEX {
    tag "${meta.id}_f${frame_idx}"
    label "process_single"

    input:
    tuple val(meta), path(complex_file), val(frame_idx)
    val(native_ligand)
    val(cofactors)

    output:
    tuple val(meta), path("${meta.id}_raw_protein.pdb"), emit: raw_protein_pdb
    tuple val(meta), path("${meta.id}_raw_native.mae"), emit: raw_native_mae
    tuple val(meta), path("${meta.id}_raw_native.sdf"), emit: raw_native_sdf

    script:
    def cof_arg = cofactors ? "--cofactors ${cofactors}" : ""
    """
    schrodinger run python3 ${projectDir}/bin/extract_complex.py \\
        --input ${complex_file} \\
        --frame_index ${frame_idx} \\
        --ligand ${native_ligand} \\
        ${cof_arg} \\
        --out_protein ${meta.id}_raw_protein.pdb \\
        --out_ligand ${meta.id}_raw_native.mae
    """

    stub:
    """
    touch ${meta.id}_raw_protein.pdb
    touch ${meta.id}_raw_native.mae
    touch ${meta.id}_raw_native.sdf
    """
}
