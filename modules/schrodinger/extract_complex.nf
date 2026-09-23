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
    val(membrane_enabled)

    output:
    tuple val(meta), path("${meta.id}_pure_protein.pdb"), emit: pure_protein_pdb, optional: true
    tuple val(meta), path("${meta.id}_complex.mae"), emit: complex_mae, optional: true
    tuple val(meta), path("${meta.id}_raw_receptor.mae"), emit: raw_receptor_mae, optional: true
    tuple val(meta), path("${meta.id}_ref_native.mae"), emit: ref_native_mae, optional: true
    tuple val(meta), path("${meta.id}_ref_native.sdf"), emit: ref_native_sdf, optional: true

    script:
    def cof_arg = cofactors ? "--cofactors ${cofactors}" : ""
    if (membrane_enabled) {
        """
        schrodinger run ${params.workdir}/bin/extract_complex.py \\
            --input ${complex_file} \\
            --frame_index ${frame_idx} \\
            --ligand ${native_ligand} \\
            ${cof_arg} \\
            --out_pure_protein ${meta.id}_pure_protein.pdb \\
            --out_complex_mae ${meta.id}_complex.mae
        """
    } else {
        """
        schrodinger run ${params.workdir}/bin/extract_complex.py \\
            --input ${complex_file} \\
            --frame_index ${frame_idx} \\
            --ligand ${native_ligand} \\
            ${cof_arg} \\
            --out_receptor_mae ${meta.id}_raw_receptor.mae \\
            --out_ligand_mae ${meta.id}_ref_native.mae
        """
    }

    stub:
    """
    touch ${meta.id}_pure_protein.pdb
    touch ${meta.id}_complex.mae
    touch ${meta.id}_raw_receptor.mae
    touch ${meta.id}_ref_native.mae
    touch ${meta.id}_ref_native.sdf
    """
}
