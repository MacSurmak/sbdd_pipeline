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
    python3 ${projectDir}/bin/wrap_gromacs.py trjcat \\
        --xtcs md_*.xtc \\
        --out_xtc concatenated.xtc
    """

    stub:
    """
    touch concatenated.xtc
    """
}

