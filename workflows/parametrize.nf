/*
 * workflows/parametrize.nf
 * On-demand chemical substance parametrization workflow supporting ORCA RESP and AM1-BCC
 */

include { PARAMETRIZE_SUBSTANCE_RESP; PARAMETRIZE_SUBSTANCE_BCC } from '../modules/amber/antechamber'

workflow PARAMETRIZE_WF {
    if (!params.param_plan) { error "FATAL: Parameter --param_plan (JSON) is required for on-demand parametrization." }

    def planFile = file(params.param_plan)
    if (!planFile.exists()) { error "FATAL: Parameter plan file not found: ${params.param_plan}" }

    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    ch_substances = Channel.fromList(planData).map { item ->
        def sdf_file = file(item.sdf)
        if (!sdf_file.exists()) error "FATAL: Reference SDF not found for ${item.name}: ${item.sdf}"

        def target_method = item.method ?: params.param_method ?: "resp"
        def meta = [
            id: item.name,
            name: item.name,
            res_name: item.res_name ?: "LIG",
            inchikey: item.inchikey,
            state_hash: item.state_hash ?: item.inchikey,
            charge: item.charge as int,
            method: target_method.toLowerCase(),
            project: params.project_name
        ]
        [ meta, sdf_file.toAbsolutePath() ]
    }

    ch_branched = ch_substances.branch {
        resp: it[0].method == 'resp'
        bcc:  it[0].method == 'bcc'
        other: true
    }

    ch_branched.other.subscribe { meta, sdf ->
        error "FATAL: Unsupported parametrization method '${meta.method}' declared for molecule ${meta.name}!"
    }

    // Branch 1: High-precision Two-Stage ORCA RESP
    PARAMETRIZE_SUBSTANCE_RESP(
        ch_branched.resp,
        params.param_ff ?: "gaff2",
        params.orca_maxcore ?: 4000,
        params.orca_opt_theory ?: "B3LYP",
        params.orca_opt_basis ?: "def2-SVP",
        params.orca_opt_dispersion ?: "D4",
        params.orca_opt_solvent ?: "CPCM(Water)",
        params.orca_opt_rijcosx ?: "true",
        params.orca_opt_maxiter ?: 200,
        params.orca_scf_conv ?: "TightSCF",
        params.orca_scf_maxiter ?: 250,
        params.orca_esp_theory ?: "HF",
        params.orca_esp_basis ?: "6-31G*",
        params.orca_esp_density ?: "Merz-Kollman"
    )

    // Branch 2: Fast Semi-empirical AM1-BCC
    PARAMETRIZE_SUBSTANCE_BCC(
        ch_branched.bcc
    )

    // Unify topology outputs into canonical 3-element channel [meta, mol2, frcmod]
    ch_resp_top = PARAMETRIZE_SUBSTANCE_RESP.out.topology.map { meta, mol2, frcmod, opt_sdf ->
        [ meta, mol2, frcmod ]
    }
    ch_all_topology = ch_resp_top.mix(PARAMETRIZE_SUBSTANCE_BCC.out.topology)

    ch_all_topology.map { meta, mol2, frcmod ->
        def state_dir = file("${params.outdir}/substances/${meta.inchikey}/states/${meta.state_hash}").toAbsolutePath()
        "${meta.state_hash}\t${state_dir}/${mol2.name}\t${state_dir}/${frcmod.name}"
    }.collectFile(name: 'substance_param_receipt.tsv', storeDir: params.outdir, newLine: true)
}
