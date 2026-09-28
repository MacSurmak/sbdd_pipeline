/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow with ensemble support
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS  } from '../subworkflows/prepare_ligands'

workflow PREPARE_WF {
    if (!params.receptor) { error "Parameter --receptor is required." }

    ch_rec = Channel.fromPath(file(params.receptor).toAbsolutePath(), checkIfExists: true)
        .map { rec -> [ [id: params.project_name, project: params.project_name], rec ] }

    PREPARE_RECEPTOR(
        ch_rec,
        params.model_index,
        params.native_ligand,
        params.keep_cofactors,
        params.membrane,
        params.membrane_type,
        params.membrane_topo,
        params.membrane_chains,
        params.target_ph,
        params.ph_tolerance,
        params.minimize_rmsd,
        params.cap_termini
    )

    // 2. Prepare Ligands (if provided)
    if (params.library) {
        ch_lib = Channel.fromPath(file(params.library).toAbsolutePath(), checkIfExists: true)
            .map { lib -> [ [id: params.project_name, project: params.project_name], lib ] }

        PREPARE_LIGANDS(
            ch_lib,
            PREPARE_RECEPTOR.out.ref_native_mae,
            params.dock_engine,
            params.batch_size,
            params.ligand_ph,
            params.ligand_ph_tol,
            params.ligprep_ff,
            params.max_states
        )
    }
}
