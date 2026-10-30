/*
 * modules/schrodinger/ligprep.nf
 * Ligand Preparation with Epik and ForceField
 */

process LIGPREP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(spec_json), path(input_sdfs, stageAs: 'in_sdfs_*/*')
    val(target_ph)
    val(ph_tolerance)
    val(ff_version)
    val(max_states)
    val(outdir)

    output:
    path("substances/**"), emit: published_substances
    path("batch_prep_receipt.tsv"), emit: receipt

    script:
    """
    schrodinger run python3 ${projectDir}/bin/wrap_schrodinger.py ligprep \\
        --spec_json ${spec_json} \\
        --in_sdfs ${input_sdfs} \\
        --outdir "${outdir}" \\
        --out_receipt batch_prep_receipt.tsv \\
        --target_ph ${target_ph} \\
        --ph_tolerance ${ph_tolerance} \\
        --ff_version ${ff_version} \\
        --max_states ${max_states}
    """
}

