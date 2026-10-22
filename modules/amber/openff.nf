/*
 * modules/amber/openff.nf
 * Builds OpenFF Sage topology (mol2 + frcmod) with guaranteed ahead-of-time disjoint namespaces
 */

process PARAMETRIZE_SUBSTANCE_OPENFF {
    tag "${meta.name}"
    label 'process_single'

    input:
    tuple val(meta), path(base_sdf), path(charges_crg)

    output:
    tuple val(meta), path("${meta.res_name}_topology.mol2"), path("${meta.res_name}_topology.frcmod"), emit: topology

    script:
    def charges_arg = (charges_crg && charges_crg.name != 'NO_CHARGES') ? "--charges ${charges_crg}" : ""
    def prefix_arg = meta.atom_type_prefix ? "--prefix ${meta.atom_type_prefix}" : "--prefix 1"
    def ff_file = meta.forcefield ? (meta.forcefield.endsWith('.offxml') ? meta.forcefield : "${meta.forcefield}.offxml") : "openff-2.3.0.offxml"
    """
    python3 ${projectDir}/bin/build_openff_topology.py \\
        --sdf ${base_sdf} \\
        --charge ${meta.charge} \\
        --res_name "${meta.res_name}" \\
        --forcefield "${ff_file}" \\
        --charge_method "${meta.method}" \\
        ${charges_arg} \\
        ${prefix_arg} \\
        --out_mol2 "${meta.res_name}_topology.mol2" \\
        --out_frcmod "${meta.res_name}_topology.frcmod"
    """
}
