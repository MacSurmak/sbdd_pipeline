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
    """
    mult=\$(python3 -c "import json; print(json.load(open('${metadata_json}'))['spin_multiplicity'])")
    python3 ${projectDir}/bin/orca_opt.py \\
        --template "${projectDir}/assets/orca/refine.inp" \\
        --input ${preopt_sdf} \\
        --charge ${meta.charge} \\
        --mult \${mult} \\
        --nprocs ${task.cpus} \\
        --maxcore_mb ${params.orca_maxcore} \\
        --out_sdf "${meta.res_name}_refined.sdf" \\
        --out_log "orca_refine.out"
    """
}
