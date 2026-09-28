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
    count=\$(ls -v md_*.xtc 2>/dev/null | wc -l)
    if [ "\$count" -eq 1 ]; then
        cp md_*.xtc concatenated.xtc
    else
        ls -v md_*.xtc > chunks.list
        gmx trjcat -f \$(cat chunks.list) -cat -o concatenated.xtc
    fi
    """

    stub:
    """
    touch concatenated.xtc
    """
}

