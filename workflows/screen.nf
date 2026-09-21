/*
 * workflows/screen.nf
 * Standalone distributed docking and screening workflow
 */

include { STAGE_INPUTS      } from "../modules/common/stage_inputs"
include { DOCK_LIBRARY      } from "../subworkflows/dock_library"
include { PREPARE_LIGANDS   } from "../subworkflows/prepare_ligands"
include { SUMMARIZE_DOCKING } from "../modules/reporting/summarize"

workflow SCREEN_WF {
    def prepared_pdb = file("${params.outdir}/receptor/${params.project_name}_prepared.pdb")
    def prepared_mae = file("${params.outdir}/receptor/${params.project_name}_prepared.mae")
    def ref_sdf      = file("${params.outdir}/receptor/${params.project_name}_ref_native.sdf")
    def ref_mae      = file("${params.outdir}/receptor/${params.project_name}_ref_native.mae")

    def rec_file = prepared_pdb.exists() ? prepared_pdb : file(params.receptor)
    def ref_file = (ref_sdf.exists() ? ref_sdf : (params.ref_native ? file(params.ref_native) : rec_file))
    def lib_file = file(params.library)

    STAGE_INPUTS(rec_file, lib_file, ref_file)

    ch_rec = STAGE_INPUTS.out.staged_receptor.map { r -> [ [id: params.project_name], r ] }
    ch_ref = STAGE_INPUTS.out.staged_ref.map { ref -> [ [id: params.project_name], ref ] }
    ch_lib = STAGE_INPUTS.out.staged_library.map { l -> [ [id: params.project_name], l ] }

    // Split library into parallel chunks
    ch_rec_mae = prepared_mae.exists() ? Channel.of([ [id: params.project_name], prepared_mae ]) : ch_rec
    ch_ref_mae = ref_mae.exists() ? Channel.of([ [id: params.project_name], ref_mae ]) : ch_ref

    // Prepare ligands (LigPrep and/or SDF chunks)
    PREPARE_LIGANDS(
        ch_lib,
        params.dock_engine ?: "gnina",
        params.batch_size ?: 15,
        7.4,
        2.0,
        16,
        32
    )

    // Parallel docking across Slurm GPU nodes
    // Parallel docking (Glide, GNINA, or Both)
    DOCK_LIBRARY(
        ch_rec,
        ch_ref,
        PREPARE_LIGANDS.out.chunks,
        ch_rec_mae,
        ch_ref_mae,
        PREPARE_LIGANDS.out.prepared_ligands,
        params.dock_engine ?: "gnina",
        params.exhaustiveness ?: 16,
        params.cnn_scoring ?: "rescore",
        params.num_modes ?: 3,
        params.autobox_add ?: 4.0
    )

    // Summarize raw poses into ranking
    SUMMARIZE_DOCKING(
        DOCK_LIBRARY.out.summary_bundle
    )
}
