/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow with ensemble support
 */

include { STAGE_INPUTS     } from '../modules/common/stage_inputs'
include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS  } from '../subworkflows/prepare_ligands'

include { COUNT_MODELS     } from '../modules/common/count_models'

workflow PREPARE_WF {
    def rec_file = params.receptor ? file(params.receptor) : file('NO_FILE')
    def lib_file = params.library ? file(params.library) : file('NO_FILE')

    // 1. Stage input files to /shared scratch
    STAGE_INPUTS(rec_file, lib_file, file('NO_FILE'))

    // 2. Dynamic determination of NMR ensemble frames
    ch_rec_meta = STAGE_INPUTS.out.staged_receptor.map { rec ->
        [ [id: params.project_name, project: params.project_name], rec ]
    }
    COUNT_MODELS(ch_rec_meta)

    ch_rec = STAGE_INPUTS.out.staged_receptor
        .combine(COUNT_MODELS.out.count_file)
        .flatMap { staged_rec, meta_cnt, count_txt ->
            def total = count_txt.text.trim() as int
            def req = params.model_index ? (params.model_index as int) : 0
            def models = (req > 0 && req <= total) ? [ req ] : (1..total).toList()
            models.collect { idx ->
                def meta = [
                    id: "${params.project_name}_f${idx}",
                    frame: idx,
                    project: params.project_name
                ]
                [ meta, staged_rec, idx ]
            }
        }

    PREPARE_RECEPTOR(
        ch_rec,
        params.native_ligand,
        params.keep_cofactors ?: "",
        params.membrane ?: false,
        params.membrane_type ?: "MOM",
        params.membrane_topo ?: "in",
        params.membrane_chains ?: "A",
        params.target_ph ?: 7.4,
        params.ph_tolerance ?: 1.0,
        params.minimize_rmsd ?: 0.3,
        params.cap_termini ?: true
    )

    // 3. Prepare Ligands (if provided)
    if (params.library) {
        ch_lib = STAGE_INPUTS.out.staged_library
            .map { staged_lib -> [ [id: params.project_name, project: params.project_name], staged_lib ] }

        PREPARE_LIGANDS(
            ch_lib,
            params.dock_engine ?: "gnina",
            params.batch_size ?: 15,
            params.ligand_ph ?: 7.4,
            params.ligand_ph_tol ?: 2.0,
            params.ligprep_ff ?: 16,
            params.max_states ?: 32
        )
    }
}
