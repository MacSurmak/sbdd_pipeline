/*
 * modules/amber/antechamber.nf
 * Strictly isolated module for AM1-BCC charge calculation and GAFF2 parametrization
 */

process PARAMETRIZE_POSE {
    tag "${meta.name}"
    label 'process_single'

    input:
    tuple val(meta), path(sdf)

    output:
    tuple val(meta), path("${meta.res_name}_pose.sdf"), path("${meta.res_name}_topology.mol2"), path("${meta.res_name}_topology.frcmod")

    script:
    """
    cp "${sdf}" "${meta.res_name}_pose.sdf"

    python3 ${projectDir}/bin/wrap_amber.py antechamber \\
        --sdf "${sdf}" \\
        --charge ${meta.charge} \\
        --mol_id "${meta.res_name}_topology" \\
        --res_name "${meta.res_name}"

    mv "${meta.res_name}_topology_gaff2.mol2" "${meta.res_name}_topology.mol2"
    """
}
