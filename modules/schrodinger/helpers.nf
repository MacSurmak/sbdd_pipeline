process CONVERT_READY_LIGAND {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(sdf, stageAs: 'in_ligand/*'), val(idx)

    output:
    tuple val(meta), path("ligand.maegz"), emit: prepared_maegz
    tuple val(meta), path("ligand.sdf"), emit: prepared_sdf

    script:
    def lig_title = meta.ligand_title ?: meta.name ?: meta.id
    """
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py convert_ligand \\
        --in_sdf "\$(ls in_ligand/* | head -n 1)" \\
        --out_maegz ligand.maegz \\
        --out_sdf ligand.sdf \\
        --title "${lig_title}"
    """
}
