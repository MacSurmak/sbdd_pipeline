/*
 * modules/membrane/packmol.nf
 * Lipid Bilayer Assembly via packmol-memgen
 */

process PACKMOL_MEMGEN {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(complex_pdb)
    val(upper_lipids)
    val(upper_ratio)
    val(lower_lipids)
    val(lower_ratio)
    val(lipid_pad)
    val(wat_pad)
    val(salt_conc)

    output:
    tuple val(meta), path("${meta.id}_mem_complex.pdb"), emit: mem_complex_pdb

    script:
    """
    python3 ${projectDir}/bin/wrap_membrane.py packmol \\
        --pdb ${complex_pdb} \\
        --upper_lipids "${upper_lipids}" \\
        --upper_ratio "${upper_ratio}" \\
        --lower_lipids "${lower_lipids}" \\
        --lower_ratio "${lower_ratio}" \\
        --lipid_pad ${lipid_pad} \\
        --wat_pad ${wat_pad} \\
        --salt_conc ${salt_conc} \\
        --out_id ${meta.id}
    """

    stub:
    """
    touch ${meta.id}_mem_complex.pdb
    """
}

