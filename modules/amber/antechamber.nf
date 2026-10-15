/*
 * modules/amber/antechamber.nf
 * Strictly isolated module for AM1-BCC charge calculation and GAFF2 parametrization
 */

process PARAMETRIZE_POSE {
    tag "${meta.name}"
    label 'process_single'
    publishDir path: { "${params.outdir}/prepared_poses/${meta.name}" }, mode: 'copy'

    input:
    tuple val(meta), path(sdf)

    output:
    tuple val(meta), path("${meta.name}_gaff2.mol2"), path("${meta.name}.frcmod")

    script:
    """
    python3 ${projectDir}/bin/wrap_amber.py antechamber \\
        --sdf "${sdf}" \\
        --charge ${meta.charge} \\
        --mol_id "${meta.name}" \\
        --res_name "${meta.res_name}"
    """

    stub:
    """
    touch "${meta.name}_gaff2.mol2"
    touch "${meta.name}.frcmod"
    """
}
