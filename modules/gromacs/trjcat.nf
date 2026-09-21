/*
 * modules/gromacs/trjcat.nf
 * Concatenates sequential MD chunk trajectories
 */

process TRJCAT {
    tag "${meta.id}_rep${meta.rep}"
    label 'process_low'

    input:
    tuple val(meta), path(xtc_chunks)

    output:
    tuple val(meta), path("concatenated.xtc"), emit: concat_xtc

    script:
    """
    # Sort files naturally by chunk index
    ls -v md_*.xtc > chunks.list
    gmx trjcat -f \$(cat chunks.list) -cat -o concatenated.xtc
    """

    stub:
    """
    touch concatenated.xtc
    """
}

