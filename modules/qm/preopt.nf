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
    python3 ${projectDir}/bin/orca_opt.py \\
        --template "${projectDir}/assets/orca/preopt.inp" \\
        --input ${valid_sdf} \\
        --charge ${meta.charge} \\
        --metadata_json ${metadata_json} \\
        --nprocs ${task.cpus} \\
        --maxcore_mb ${params.orca_maxcore} \\
        --out_sdf "${meta.res_name}_preopt.sdf" \\
        --out_log "orca_preopt.out"
    """
}
