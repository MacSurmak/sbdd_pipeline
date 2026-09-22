/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow with ensemble support
 */

include { STAGE_INPUTS     } from '../modules/common/stage_inputs'
include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS  } from '../subworkflows/prepare_ligands'

workflow PREPARE_WF {
    def rec_file = params.receptor ? file(params.receptor) : file('NO_FILE')
    def lib_file = params.library ? file(params.library) : file('NO_FILE')

    // 1. Stage input files from local home to /shared scratch
    STAGE_INPUTS(rec_file, lib_file, file('NO_FILE'))

    // 2. Determine models channel:
    // If params.model_index > 0: run single specified model.
    // If params.model_index == 0 or null: run all ensemble models.
    def requested_model = params.model_index ? (params.model_index as int) : 0
    def num_models = requested_model > 0 ? requested_model : (params.total_models ? (params.total_models as int) : 20)
    def model_list = requested_model > 0 ? [ requested_model ] : (1..num_models).toList()

    ch_rec = STAGE_INPUTS.out.staged_receptor
        .combine(Channel.fromList(model_list))
        .map { staged_rec, idx ->
            def meta = [
                id: "${params.project_name}_f${idx}",
                frame: idx,
                project: params.project_name
            ]
            [ meta, staged_rec, idx ]
        }

    PREPARE_RECEPTOR(
        ch_rec,
        params.native_ligand,
        params.keep_cofactors ?: "",
        params.membrane ?: false,
        params.membrane_type ?: "MOM",
        params.membrane_topo ?: "in",
        params.membrane_chains ?: "A",
        7.4,
        1.0,
        0.3,
        true
    )

    // 3. Prepare Ligands (if provided)
    if (params.library) {
        ch_lib = STAGE_INPUTS.out.staged_library
            .map { staged_lib -> [ [id: params.project_name, project: params.project_name], staged_lib ] }

        PREPARE_LIGANDS(
            ch_lib,
            params.dock_engine ?: "gnina",
            params.batch_size ?: 15,
            7.4,
            2.0,
            16,
            32
        )
    }
}
