/*
 * workflows/screen.nf
 * Standalone distributed docking and screening workflow with ensemble support
 */

include { STAGE_INPUTS      } from "../modules/common/stage_inputs"
include { PREPARE_RECEPTOR  } from "../subworkflows/prepare_receptor"
include { PREPARE_LIGANDS   } from "../subworkflows/prepare_ligands"
include { DOCK_LIBRARY      } from "../subworkflows/dock_library"
include { SUMMARIZE_DOCKING } from "../modules/reporting/summarize"

workflow SCREEN_WF {
    def rec_file = params.receptor ? file(params.receptor) : file('NO_FILE')
    def lib_file = params.library ? file(params.library) : file('NO_FILE')

    // Discover or generate all ensemble frames
    def requested_model = params.model_index ? (params.model_index as int) : 0
    def num_models = requested_model > 0 ? requested_model : (params.total_models ? (params.total_models as int) : 20)
    def model_list = requested_model > 0 ? [ requested_model ] : (1..num_models).toList()

    def prep_check = file("${params.outdir}/receptor/${params.project_name}_f1_prepared.pdb")

    if (prep_check.exists()) {
        ch_rec_pdb = Channel.fromList(model_list).map { idx ->
            def meta = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
            [ meta, file("${params.outdir}/receptor/${params.project_name}_f${idx}_prepared.pdb") ]
        }
        ch_ref_sdf = Channel.fromList(model_list).map { idx ->
            def meta = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
            [ meta, file("${params.outdir}/receptor/${params.project_name}_f${idx}_ref_native.sdf") ]
        }
        ch_rec_mae = Channel.fromList(model_list).map { idx ->
            def meta = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
            [ meta, file("${params.outdir}/receptor/${params.project_name}_f${idx}_prepared.mae") ]
        }
        ch_ref_mae = Channel.fromList(model_list).map { idx ->
            def meta = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
            [ meta, file("${params.outdir}/receptor/${params.project_name}_f${idx}_ref_native.mae") ]
        }

        STAGE_INPUTS(rec_file, lib_file, file('NO_FILE'))
    } else {
        STAGE_INPUTS(rec_file, lib_file, file('NO_FILE'))

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

        ch_rec_pdb = PREPARE_RECEPTOR.out.prepared_pdb
        ch_ref_sdf = PREPARE_RECEPTOR.out.ref_native_sdf
        ch_rec_mae = PREPARE_RECEPTOR.out.prepared_mae
        ch_ref_mae = PREPARE_RECEPTOR.out.ref_native_mae
    }

    // Prepare ligands (LigPrep and/or SDF chunks)
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

    // Parallel docking (Glide, GNINA, or Both)
    DOCK_LIBRARY(
        ch_rec_pdb,
        ch_ref_sdf,
        PREPARE_LIGANDS.out.chunks,
        ch_rec_mae,
        ch_ref_mae,
        PREPARE_LIGANDS.out.prepared_ligands,
        params.dock_engine ?: "gnina",
        params.exhaustiveness ?: 16,
        params.cnn_scoring ?: "rescore",
        params.num_modes ?: 3,
        params.autobox_add ?: 4.0,
        params.top_k_frames ?: 3,
        params.percentile_cutoff ?: 0.50,
        params.membrane ?: false
    )

    // Summarize raw poses into consensus ranking CSV tables
    SUMMARIZE_DOCKING(
        params.project_name,
        DOCK_LIBRARY.out.summary_files
    )
}
