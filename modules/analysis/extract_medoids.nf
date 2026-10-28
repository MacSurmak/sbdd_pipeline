/*
 * modules/analysis/extract_medoids.nf
 * Distributed Slurm process: extracts medoid complexes (Receptor, Ligand, Cofactors)
 * directly in $SLURM_TMPDIR without master node I/O bottlenecks.
 */

process EXTRACT_MEDOID_COMPLEXES {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(ref_pdb), path(xtc_files, stageAs: 'trajs_*/*'), path(spec_json)

    output:
    path("targets/**"), emit: published_frames
    path("substances/**"), emit: published_substances
    path("extract_medoids_receipt.tsv"), emit: receipt

    script:
    """
    python3 ${projectDir}/bin/extract_medoids.py \\
        --spec ${spec_json} \\
        --outdir .
    """
}
