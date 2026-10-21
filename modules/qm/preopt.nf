/*
 * modules/qm/preopt.nf
 * Stage 1A: Fast pure GGA geometry pre-relaxation (PBE-D4/def2-SVP/CPCM) without exact exchange
 */

process ORCA_PREOPT {
    tag "${meta.name}"

    input:
    tuple val(meta), path(valid_sdf), path(metadata_json)

    output:
    tuple val(meta), path("${meta.res_name}_preopt.sdf"), path(metadata_json), emit: preopt_sdf

    script:
    """
    mult=\$(python3 -c "import json; print(json.load(open('${metadata_json}'))['spin_multiplicity'])")
    python3 ${projectDir}/bin/orca_opt.py \\
        --input ${valid_sdf} \\
        --charge ${meta.charge} \\
        --mult \${mult} \\
        --theory "${params.orca_preopt_theory}" \\
        --basis "${params.orca_preopt_basis}" \\
        --dispersion "${params.orca_opt_dispersion}" \\
        --solvent "${params.orca_opt_solvent}" \\
        --opt_level "${params.orca_preopt_opt_level}" \\
        --maxiter ${params.orca_opt_maxiter} \\
        --scf_convergence "${params.orca_geom_scf_conv}" \\
        --scf_maxiter ${params.orca_geom_scf_maxiter} \\
        --nprocs ${task.cpus} \\
        --maxcore_mb ${params.orca_maxcore} \\
        --out_sdf "${meta.res_name}_preopt.sdf" \\
        --out_log "orca_preopt.out"
    """
}
