/*
 * modules/amber/build_topology.nf
 * Builds GAFF2 mol2 and frcmod parameter files using pre-calculated atomic charges
 */

process BUILD_TOPOLOGY_GAFF2 {
    tag "${meta.name}"
    label 'process_single'

    input:
    tuple val(meta), path(base_sdf), path(charges_crg)

    output:
    tuple val(meta), path("${meta.res_name}_topology.mol2"), path("${meta.res_name}_topology.frcmod"), emit: topology

    script:
    """
    python3 ${projectDir}/bin/build_gaff2_topology.py \\
        --sdf ${base_sdf} \\
        --charges ${charges_crg} \\
        --charge ${meta.charge} \\
        --res_name "${meta.res_name}" \\
        --forcefield "${params.param_ff}" \\
        --out_mol2 "${meta.res_name}_topology.mol2" \\
        --out_frcmod "${meta.res_name}_topology.frcmod"
    """
}
