/*
 * modules/amber/antechamber.nf
 * Strictly isolated module for AM1-BCC charge calculation and GAFF2 parametrization
 */

process PARAMETRIZE_POSE {
    tag "${meta.name}"
    label 'process_single'

    input:
    tuple val(meta), path(sdf), path(mae)

    output:
    tuple val(meta), path("${meta.res_name}_pose.sdf"), path("${meta.res_name}_pose.mae", optional: true), path("${meta.res_name}_topology.mol2", optional: true), path("${meta.res_name}_topology.frcmod", optional: true)

    script:
    if (meta.is_ion) {
        """
        cp "${sdf}" "${meta.res_name}_pose.sdf"
        if [ -f "${mae}" ] && [ -s "${mae}" ]; then cp "${mae}" "${meta.res_name}_pose.mae"; fi
        """
    } else {
        """
        cp "${sdf}" "${meta.res_name}_pose.sdf"
        if [ -f "${mae}" ] && [ -s "${mae}" ]; then cp "${mae}" "${meta.res_name}_pose.mae"; fi

        python3 ${projectDir}/bin/wrap_amber.py antechamber \\
            --sdf "${sdf}" \\
            --charge ${meta.charge} \\
            --mol_id "${meta.res_name}_topology" \\
            --res_name "${meta.res_name}"

        mv "${meta.res_name}_topology_gaff2.mol2" "${meta.res_name}_topology.mol2"
        """
    }
}
