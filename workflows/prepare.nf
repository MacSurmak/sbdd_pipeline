/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow
 */

include { STAGE_INPUTS     } from '../modules/common/stage_inputs'
include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS  } from '../subworkflows/prepare_ligands'

workflow PREPARE_WF {
    // Receptor
    ch_rec = Channel.of([ [id: params.project_name], file(params.receptor), params.model_index ?: 1 ])
    
    def rec_file = file(params.receptor)
    def lib_file = params.library ? file(params.library) : file('NO_FILE')

    // 1. Stage input files from local home to /shared scratch
    STAGE_INPUTS(rec_file, lib_file, file('NO_FILE'))

    // 2. Prepare Receptor
    ch_rec = STAGE_INPUTS.out.staged_receptor
        .map { staged_rec -> [ [id: params.project_name], staged_rec, params.model_index ?: 1 ] }

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

    // Ligands (if provided)
    // 3. Prepare Ligands (if provided)
    if (params.library) {
        ch_lib = Channel.of([ [id: params.project_name], file(params.library) ])
        ch_lib = STAGE_INPUTS.out.staged_library
            .map { staged_lib -> [ [id: params.project_name], staged_lib ] }

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

