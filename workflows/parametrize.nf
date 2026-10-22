/*
 * workflows/parametrize.nf
 * On-demand chemical substance parametrization workflow supporting ORCA RESP and AM1-BCC
 */

include { VALIDATE_LIGAND               } from '../modules/qm/validate'
include { ORCA_PREOPT                   } from '../modules/qm/preopt'
include { ORCA_REFINE_OPT               } from '../modules/qm/refine_opt'
include { ORCA_ESP_CALC                 } from '../modules/qm/esp'
include { FIT_RESP_CHARGES              } from '../modules/qm/charges'
include { BUILD_TOPOLOGY_GAFF2          } from '../modules/amber/build_topology'
include { PARAMETRIZE_SUBSTANCE_BCC     } from '../modules/amber/antechamber'
include { PARAMETRIZE_SUBSTANCE_OPENFF  } from '../modules/amber/openff'

workflow PARAMETRIZE_WF {
    if (!params.param_plan) { error "FATAL: Parameter --param_plan (JSON) is required for on-demand parametrization." }

    def planFile = file(params.param_plan)
    if (!planFile.exists()) { error "FATAL: Parameter plan file not found: ${params.param_plan}" }

    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    ch_substances = Channel.fromList(planData).map { item ->
        def sdf_file = file(item.sdf)
        if (!sdf_file.exists()) error "FATAL: Reference SDF not found for ${item.name}: ${item.sdf}"

        def target_method = item.method ?: params.param_method ?: "ashgc"
        def target_ff = item.forcefield ?: params.param_ff ?: "openff-2.3.0"
        def meta = [
            id: item.name,
            name: item.name,
            res_name: item.res_name ?: "LIG",
            inchikey: item.inchikey,
            state_hash: item.state_hash ?: item.inchikey,
            charge: item.charge as int,
            method: target_method.toLowerCase(),
            forcefield: target_ff.toLowerCase(),
            atom_type_prefix: item.atom_type_prefix ?: "1",
            project: params.project_name
        ]
        [ meta, sdf_file.toAbsolutePath() ]
    }

    def method = params.param_method ? params.param_method.toLowerCase() : "ashgc"
    def ff = params.param_ff ? params.param_ff.toLowerCase() : "openff-2.3.0"

    if (ff.startsWith("openff")) {
        if (method == "resp") {
            VALIDATE_LIGAND(ch_substances)
            ORCA_PREOPT(VALIDATE_LIGAND.out.validated)
            ORCA_REFINE_OPT(ORCA_PREOPT.out.preopt_sdf)
            ORCA_ESP_CALC(ORCA_REFINE_OPT.out.refined_sdf)
            FIT_RESP_CHARGES(ORCA_ESP_CALC.out.esp_out)

            ch_openff_in = ch_substances.join(FIT_RESP_CHARGES.out.charges, by: 0)
            PARAMETRIZE_SUBSTANCE_OPENFF(ch_openff_in)
            ch_topology = PARAMETRIZE_SUBSTANCE_OPENFF.out.topology
        } else if (method in ["ashgc", "bcc"]) {
            ch_dummy_charges = file("${workDir}/NO_CHARGES")
            ch_dummy_charges.text = ""
            ch_openff_in = ch_substances.map { meta, sdf -> [ meta, sdf, ch_dummy_charges ] }
            PARAMETRIZE_SUBSTANCE_OPENFF(ch_openff_in)
            ch_topology = PARAMETRIZE_SUBSTANCE_OPENFF.out.topology
        } else {
            error "FATAL: Unsupported parametrization method '${method}' for OpenFF! Must be 'ashgc', 'resp', or 'bcc'."
        }
    } else if (ff in ["gaff", "gaff2"]) {
        if (method == "resp") {
            VALIDATE_LIGAND(ch_substances)
            ORCA_PREOPT(VALIDATE_LIGAND.out.validated)
            ORCA_REFINE_OPT(ORCA_PREOPT.out.preopt_sdf)
            ORCA_ESP_CALC(ORCA_REFINE_OPT.out.refined_sdf)
            FIT_RESP_CHARGES(ORCA_ESP_CALC.out.esp_out)

            ch_build_in = ch_substances.join(FIT_RESP_CHARGES.out.charges, by: 0)
            BUILD_TOPOLOGY_GAFF2(ch_build_in)
            ch_topology = BUILD_TOPOLOGY_GAFF2.out.topology
        } else if (method == "bcc") {
            PARAMETRIZE_SUBSTANCE_BCC(ch_substances)
            ch_topology = PARAMETRIZE_SUBSTANCE_BCC.out.topology
        } else {
            error "FATAL: Unsupported parametrization method '${method}' for GAFF2! Must be 'resp' or 'bcc'."
        }
    } else {
        error "FATAL: Unsupported force field '${ff}' declared in configuration! Must be 'openff-2.3.0', 'openff-2.2.1', 'gaff2', or 'gaff'."
    }

    ch_topology.map { meta, mol2, frcmod ->
        def state_dir = file("${params.outdir}/substances/${meta.inchikey}/states/${meta.state_hash}").toAbsolutePath()
        "${meta.state_hash}\t${state_dir}/${mol2.name}\t${state_dir}/${frcmod.name}"
    }.collectFile(name: 'substance_param_receipt.tsv', storeDir: params.outdir, newLine: true)
}
