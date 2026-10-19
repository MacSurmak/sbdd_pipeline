/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow with ensemble support
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS  } from '../subworkflows/prepare_ligands'

workflow PREPARE_WF {
    if (!params.plan) { error "FATAL: Parameter --plan (JSON) is required for database-driven preparation." }
    
    def planFile = file(params.plan)
    if (!planFile.exists()) { error "FATAL: Plan file not found at ${params.plan}" }
    
    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    // 1. Process Receptors if present in plan
    def receptorList = planData.receptors ?: (planData instanceof List ? planData : [])

    if (receptorList) {
        ch_targets = Channel.fromList(receptorList).map { target ->
            def rec_path = file(target.path)
            if (!rec_path.exists()) {
                error "FATAL: Path explicitly provided by DB does not exist: ${target.path}"
            }

            def meta = [
                id: target.name,
                name: target.name,
                target: target.target_name,
                prep_mode: target.prep_mode,
                project: params.project_name,
                native_ligand: target.native_ligand ?: "",
                cofactors: target.cofactors ?: [],
                heterogens: target.heterogens ?: []
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

    // 2. Process Ligands if present in plan
    def ligandList = planData.ligands ?: []
    if (ligandList) {
        ch_lig_in = Channel.fromList(ligandList).map { lig ->
            def lig_path = file(lig.path)
            if (!lig_path.exists()) {
                error "FATAL: Ligand file path does not exist: ${lig.path}"
            }
            def meta = [
                id: lig.name,
                name: lig.name,
                inchikey: lig.inchikey,
                project: params.project_name,
                status: lig.status ?: 'unprepared',
                is_native: lig.is_native ?: false
            ]
            [ meta, lig_path.toAbsolutePath() ]
        }

        if (!params.batch_size || !params.ligand_ph || !params.ligand_ph_tol || !params.ligprep_ff || !params.max_states) {
            error "FATAL: Ligand preparation parameters missing in configuration! Check project.yaml."
        }

        PREPARE_LIGANDS(
            ch_lig_in,
            params.batch_size,
            params.ligand_ph,
            params.ligand_ph_tol,
            params.ligprep_ff,
            params.max_states
        )

        PREPARE_LIGANDS.out.prepared_ligands.map { meta, mae ->
            def base_id = meta.name ?: meta.id
            def sub_dir = file("${params.outdir}/substances/${meta.inchikey}").toAbsolutePath()
            "${base_id}\t${sub_dir}/${meta.name}_conformers.maegz\t${sub_dir}/${meta.name}_conformers.sdf"
        }.collectFile(name: 'prep_ligand_receipt.tsv', storeDir: params.outdir, newLine: true)
    }
}
