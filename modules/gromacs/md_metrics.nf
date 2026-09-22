/*
 * modules/gromacs/md_metrics.nf
 * Calculates trajectory geometry metrics (RMSD, RMSF, Contacts)
 * Self-contained worker process: zero external script file staging.
 */

process MD_METRICS {
    tag "${meta.id}_rep${meta.rep}"
    label "process_medium"

    input:
    tuple val(meta), path(clean_xtc), path(ref_pdb)

    output:
    tuple val(meta), path("rmsd.csv"), path("rmsf.csv"), path("contacts.csv"), emit: tables
    tuple val(meta), path("md_metrics_summary.json"), emit: summary_json

    script:
    """
    python3 ${projectDir}/bin/md_metrics.py \\
        --pdb ${ref_pdb} \\
        --xtc ${clean_xtc} \\
        --outdir .
    """

    stub:
    """
    touch rmsd.csv
    touch rmsf.csv
    touch contacts.csv
    touch md_metrics_summary.json
    """
}
