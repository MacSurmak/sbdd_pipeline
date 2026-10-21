/*
 * modules/qm/refine_opt.nf
 * Stage 1B: High-level hybrid DFT geometry refinement (B3LYP-D4/def2-SVP/RIJCOSX/CPCM)
 */

process ORCA_REFINE_OPT {
    tag "${meta.name}"

    input:
    tuple val(meta), path(preopt_sdf), path(metadata_json)

    output:
    tuple val(meta), path("${meta.res_name}_refined.sdf"), path(metadata_json), emit: refined_sdf

    script:
    def rij_flag = (params.orca_opt_rijcosx.toString().toLowerCase() in ['true', '1', 'yes']) ? "--rijcosx" : ""
    """
    mult=\$(python3 -c "import json; print(json.load(open('${metadata_json}'))['spin_multiplicity'])")
    python3 ${projectDir}/bin/orca_opt.py \\
        --input ${preopt_sdf} \\
        --charge ${meta.charge} \\
        --mult \${mult} \\
        --theory "${params.orca_refine_theory}" \\
        --basis "${params.orca_refine_basis}" \\
        --dispersion "${params.orca_opt_dispersion}" \\
        --solvent "${params.orca_opt_solvent}" \\
        ${rij_flag} \\
        --maxiter ${params.orca_opt_maxiter} \\
        --scf_convergence "${params.orca_geom_scf_conv}" \\
        --scf_maxiter ${params.orca_geom_scf_maxiter} \\
        --nprocs ${task.cpus} \\
        --maxcore_mb ${params.orca_maxcore} \\
        --out_sdf "${meta.res_name}_refined.sdf" \\
        --out_log "orca_refine.out"
    """
}
