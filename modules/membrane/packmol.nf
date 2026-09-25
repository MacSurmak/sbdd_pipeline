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
    packmol-memgen \
        --pdb ${complex_pdb} \
        --lipids "${upper_lipids}//${lower_lipids}" \
        --ratio "${upper_ratio}//${lower_ratio}" \
        --preoriented \
        --dist ${lipid_pad} \
        --dist_wat ${wat_pad} \
        --salt \
        --salt_c "K+" \
        --saltcon ${salt_conc} \
        --notprotonate \
        -o ${meta.id}_mem_complex

    # Packmol-memgen outputs <output_prefix>.pdb
    if [ -f "${meta.id}_mem_complex.pdb" ]; then
        echo "Packmol bilayer generated successfully."
    elif [ -f "${meta.id}_mem_complex_mempack.pdb" ]; then
        mv "${meta.id}_mem_complex_mempack.pdb" "${meta.id}_mem_complex.pdb"
    fi
    """

    stub:
    """
    touch ${meta.id}_mem_complex.pdb
    """
}

