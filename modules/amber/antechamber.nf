/*
 * modules/amber/antechamber.nf
 * Strictly isolated processes for small molecule parametrization (ORCA RESP and SQM AM1-BCC)
 */

process PARAMETRIZE_SUBSTANCE_RESP {
    tag "${meta.name}"
    label 'process_high'

    input:
    tuple val(meta), path(base_sdf, stageAs: 'input_base.sdf')
    val(forcefield)
    val(maxcore)
    val(opt_theory)
    val(opt_basis)
    val(opt_dispersion)
    val(opt_solvent)
    val(opt_rijcosx)
    val(opt_maxiter)
    val(scf_conv)
    val(scf_maxiter)
    val(esp_theory)
    val(esp_basis)
    val(esp_density)

    output:
    tuple val(meta), path("${meta.res_name}_topology.mol2"), path("${meta.res_name}_topology.frcmod"), path("${meta.res_name}_opt.sdf"), emit: topology

    script:
    """
    python3 ${projectDir}/bin/wrap_orca_resp.py \\
        --sdf input_base.sdf \\
        --charge ${meta.charge} \\
        --mol_id "${meta.res_name}_topology" \\
        --res_name "${meta.res_name}" \\
        --out_mol2 "${meta.res_name}_topology.mol2" \\
        --out_frcmod "${meta.res_name}_topology.frcmod" \\
        --out_opt_sdf "${meta.res_name}_opt.sdf" \\
        --forcefield "${forcefield}" \\
        --nprocs ${task.cpus} \\
        --maxcore_mb ${maxcore} \\
        --opt_theory "${opt_theory}" \\
        --opt_basis "${opt_basis}" \\
        --opt_dispersion "${opt_dispersion}" \\
        --opt_solvent "${opt_solvent}" \\
        --opt_rijcosx "${opt_rijcosx}" \\
        --opt_maxiter ${opt_maxiter} \\
        --scf_convergence "${scf_conv}" \\
        --scf_maxiter ${scf_maxiter} \\
        --esp_theory "${esp_theory}" \\
        --esp_basis "${esp_basis}" \\
        --esp_density "${esp_density}"
    """
}

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
