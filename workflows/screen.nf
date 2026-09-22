/*
 * workflows/screen.nf
 * Standalone distributed docking and screening workflow with ensemble support
 */

include { STAGE_INPUTS      } from "../modules/common/stage_inputs"
include { PREPARE_RECEPTOR  } from "../subworkflows/prepare_receptor"
include { PREPARE_LIGANDS   } from "../subworkflows/prepare_ligands"
include { DOCK_LIBRARY      } from "../subworkflows/dock_library"
include { SUMMARIZE_DOCKING } from "../modules/reporting/summarize"

include { COUNT_MODELS      } from "../modules/common/count_models"

workflow SCREEN_WF {
    def rec_file = params.receptor ? file(params.receptor) : file('NO_FILE')
    def lib_file = params.library ? file(params.library) : file('NO_FILE')

    STAGE_INPUTS(rec_file, lib_file, file('NO_FILE'))

    // Dynamic determination of NMR ensemble frames
    ch_rec_meta = STAGE_INPUTS.out.staged_receptor.map { rec ->
        [ [id: params.project_name, project: params.project_name], rec ]
    }
    COUNT_MODELS(ch_rec_meta)

    def prep_check = file("${params.outdir}/receptor/${params.project_name}_f1_prepared.pdb")

    if (prep_check.exists()) {
        ch_rec_pdb = COUNT_MODELS.out.count_file.flatMap { meta, count_txt ->
            def total = count_txt.text.trim() as int
            def req = params.model_index ? (params.model_index as int) : 0
            def models = (req > 0 && req <= total) ? [ req ] : (1..total).toList()
            models.collect { idx ->
                def m = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
                [ m, file("${params.outdir}/receptor/${params.project_name}_f${idx}_prepared.pdb") ]
            }
        }
        ch_ref_sdf = COUNT_MODELS.out.count_file.flatMap { meta, count_txt ->
            def total = count_txt.text.trim() as int
            def req = params.model_index ? (params.model_index as int) : 0
            def models = (req > 0 && req <= total) ? [ req ] : (1..total).toList()
            models.collect { idx ->
                def m = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
                [ m, file("${params.outdir}/receptor/${params.project_name}_f${idx}_ref_native.sdf") ]
            }
        }
        ch_rec_mae = COUNT_MODELS.out.count_file.flatMap { meta, count_txt ->
            def total = count_txt.text.trim() as int
            def req = params.model_index ? (params.model_index as int) : 0
            def models = (req > 0 && req <= total) ? [ req ] : (1..total).toList()
            models.collect { idx ->
                def m = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
                [ m, file("${params.outdir}/receptor/${params.project_name}_f${idx}_prepared.mae") ]
            }
        }
        ch_ref_mae = COUNT_MODELS.out.count_file.flatMap { meta, count_txt ->
            def total = count_txt.text.trim() as int
            def req = params.model_index ? (params.model_index as int) : 0
            def models = (req > 0 && req <= total) ? [ req ] : (1..total).toList()
            models.collect { idx ->
                def m = [id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name]
                [ m, file("${params.outdir}/receptor/${params.project_name}_f${idx}_ref_native.mae") ]
            }
        }
    } else {
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

        ch_rec_pdb = PREPARE_RECEPTOR.out.prepared_pdb
        ch_ref_sdf = PREPARE_RECEPTOR.out.ref_native_sdf
        ch_rec_mae = PREPARE_RECEPTOR.out.prepared_mae
        ch_ref_mae = PREPARE_RECEPTOR.out.ref_native_mae
    }

    // Prepare ligands dynamically
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
