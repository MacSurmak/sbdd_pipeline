/*
 * modules/amber/antechamber.nf
 * Strictly isolated processes for small molecule parametrization (ORCA RESP and SQM AM1-BCC)
 */

process PARAMETRIZE_SUBSTANCE_BCC {
    tag "${meta.name}"
    label 'process_single'

    input:
    tuple val(meta), path(base_sdf, stageAs: 'input_base.sdf')

    output:
    tuple val(meta), path("${meta.res_name}_topology.mol2"), path("${meta.res_name}_topology.frcmod"), emit: topology

    script:
    """
    python3 ${projectDir}/bin/wrap_amber.py antechamber \\
        --sdf input_base.sdf \\
        --charge ${meta.charge} \\
        --mol_id "${meta.res_name}_topology" \\
        --res_name "${meta.res_name}"

    mv "${meta.res_name}_topology_gaff2.mol2" "${meta.res_name}_topology.mol2"
    """
}
