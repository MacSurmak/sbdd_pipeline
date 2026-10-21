/*
 * modules/qm/validate.nf
 * Pre-flight validation of ligand 3D geometry and spin multiplicity calculation
 */

process VALIDATE_LIGAND {
    tag "${meta.name}"
    label 'process_single'

    input:
    tuple val(meta), path(base_sdf, stageAs: 'input_base.sdf')

    output:
    tuple val(meta), path("${meta.res_name}_valid.sdf"), path("metadata.json"), emit: validated

    script:
    """
    python3 ${projectDir}/bin/validate_ligand.py \\
        --sdf input_base.sdf \\
        --mol_id "${meta.res_name}" \\
        --charge ${meta.charge} \\
        --out_sdf "${meta.res_name}_valid.sdf" \\
        --out_json metadata.json
    """
}
