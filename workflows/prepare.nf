/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow with ensemble support
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'

workflow PREPARE_WF {
    if (!params.plan) { error "FATAL: Parameter --plan (JSON) is required for database-driven preparation." }
    
    def planFile = file(params.plan)
    if (!planFile.exists()) { error "FATAL: Plan file not found at ${params.plan}" }
    
    def jsonSlurper = new groovy.json.JsonSlurper()
    def planData = jsonSlurper.parseText(planFile.text)

    ch_targets = Channel.fromList(planData).map { target ->
        def rec_path = file(target.path)
        if (!rec_path.exists()) {
            error "FATAL: Path explicitly provided by DB does not exist: ${target.path}"
        }

        def meta = [
            id: target.name,
            name: target.name,
            prep_mode: target.prep_mode,
            project: params.project_name,
            native_ligand: target.native_ligand ?: "",
            cofactors: target.cofactors ?: []
        ]
        
        [ meta, rec_path.toAbsolutePath() ]
    }

    PREPARE_RECEPTOR(
        ch_targets,
        params.membrane,
        params.membrane_type,
        params.membrane_topo,
        params.membrane_chains,
        params.target_ph,
        params.ph_tolerance,
        params.minimize_rmsd,
        params.cap_termini
    )
}
