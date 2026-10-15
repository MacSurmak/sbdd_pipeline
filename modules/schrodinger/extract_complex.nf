/*
 * modules/schrodinger/extract_complex.nf
 * Extracts raw protein PDB and native complex MAE from atomic single-frame PDB/CIF.
 */

process EXTRACT_COMPLEX {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), path(complex_file)
    val(membrane_enabled)

    output:
    tuple val(meta), path("receptor_pure.pdb"), emit: pure_protein_pdb, optional: true
    tuple val(meta), path("complex_raw.mae"), emit: complex_mae

    script:
    def pure_prot_arg = (membrane_enabled.toString() == 'true') ? "--out_pure_protein receptor_pure.pdb" : ""
    """
    schrodinger run ${projectDir}/bin/extract_complex.py \\
        --input ${complex_file} \\
        --frame_index 1 \\
        ${pure_prot_arg} \\
        --out_complex_mae complex_raw.mae
    """
}
