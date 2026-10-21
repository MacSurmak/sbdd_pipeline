/*
 * modules/qm/charges.nf
 * Extraction and topological graph automorphism symmetrization of Merz-Kollman ESP charges
 */

process FIT_RESP_CHARGES {
    tag "${meta.name}"
    label 'process_single'

    input:
    tuple val(meta), path(refined_sdf), path(orca_esp_out)

    output:
    tuple val(meta), path(refined_sdf), path("${meta.res_name}_charges.crg"), emit: charges

    script:
    """
    python3 ${projectDir}/bin/fit_resp_charges.py \\
        --sdf ${refined_sdf} \\
        --orca_log ${orca_esp_out} \\
        --charge ${meta.charge} \\
        --out_crg "${meta.res_name}_charges.crg" \\
        --out_sdf "${meta.res_name}_charged.sdf"
    """
}
