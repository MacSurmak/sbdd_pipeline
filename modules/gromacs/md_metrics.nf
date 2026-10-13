/*
 * modules/gromacs/md_metrics.nf
 * Calculates trajectory geometry metrics (RMSD, RMSF, Contacts)
 * Self-contained worker process: zero external script file staging.
 */

process MD_METRICS {
    tag "${meta.id}"
    label "process_medium"

    input:
    tuple val(meta), path(clean_xtc), path(ref_pdb)
    val(c_heavy)
    val(c_hyd)
    val(c_pol)

    output:
    tuple val(meta), path("rmsd.csv"), path("rmsf.csv"), path("contacts.csv"), emit: tables
    tuple val(meta), path("md_metrics_summary.json"), emit: summary_json
    tuple val(meta), path("*.png"), emit: plots, optional: true

    script:
    def lig_arg = meta.ligand ? "--ligand ${meta.ligand}" : ""
    def cof_list = meta.cofactors ? (meta.cofactors instanceof List ? meta.cofactors.join(' ') : meta.cofactors) : ""
    def cof_arg = cof_list ? "--cofactors ${cof_list}" : ""
    """
    python3 ${projectDir}/bin/md_metrics.py \\
        --struct ${ref_pdb} \\
        --xtc ${clean_xtc} \\
        ${lig_arg} \\
        ${cof_arg} \\
        --heavy_cutoff ${c_heavy} \\
        --hydrophobic_cutoff ${c_hyd} \\
        --polar_cutoff ${c_pol} \\
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
