/*
 * workflows/parametrize.nf
 * On-demand quantum parametrization of chemical substances (GAFF2 / AM1-BCC)
 */

include { PARAMETRIZE_SUBSTANCE } from '../modules/amber/antechamber'

workflow PARAMETRIZE_WF {
    if (!params.param_plan) { error "FATAL: Parameter --param_plan (JSON) is required for on-demand parametrization." }

    def planFile = file(params.param_plan)
    if (!planFile.exists()) { error "FATAL: Parameter plan file not found: ${params.param_plan}" }

    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    ch_substances = Channel.fromList(planData).map { item ->
        def sdf_file = file(item.sdf)
        if (!sdf_file.exists()) error "FATAL: Reference SDF not found for ${item.name}: ${item.sdf}"

        def meta = [
            id: item.name,
            name: item.name,
            res_name: item.res_name ?: "LIG",
            inchikey: item.inchikey,
            state_hash: item.state_hash ?: item.inchikey,
            charge: item.charge as int,
            project: params.project_name
        ]
        [ meta, sdf_file.toAbsolutePath() ]
    }

    PARAMETRIZE_SUBSTANCE(ch_substances)

    PARAMETRIZE_SUBSTANCE.out.topology.map { meta, mol2, frcmod ->
        def state_dir = file("${params.outdir}/substances/${meta.inchikey}/states/${meta.state_hash}").toAbsolutePath()
        "${meta.state_hash}\t${state_dir}/${mol2.name}\t${state_dir}/${frcmod.name}"
    }.collectFile(name: 'substance_param_receipt.tsv', storeDir: params.outdir, newLine: true)
}
