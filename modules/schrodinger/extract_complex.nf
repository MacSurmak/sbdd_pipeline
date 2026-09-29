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
    val(membrane_enabled)

    output:
    tuple val(meta), path("${meta.id}_pure_protein.pdb"), emit: pure_protein_pdb, optional: true
    tuple val(meta), path("${meta.id}_complex.mae"), emit: complex_mae

    script:
    def pure_prot_arg = (membrane_enabled.toString() == 'true') ? "--out_pure_protein ${meta.id}_pure_protein.pdb" : ""
    """
    schrodinger run ${projectDir}/bin/extract_complex.py \\
        --input ${complex_file} \\
        --frame_index ${frame_idx} \\
        \${pure_prot_arg} \\
        --out_complex_mae ${meta.id}_complex.mae
    """
    stub:
    """
    touch ${meta.id}_pure_protein.pdb
    touch ${meta.id}_complex.mae
    """
}
