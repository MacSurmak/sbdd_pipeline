/*
 * workflows/screen.nf
 * Standalone distributed docking and screening workflow with ensemble support
 */

include { PREPARE_RECEPTOR  } from "../subworkflows/prepare_receptor"
include { PREPARE_LIGANDS   } from "../subworkflows/prepare_ligands"
include { DOCK_LIBRARY      } from "../subworkflows/dock_library"
include { SUMMARIZE_DOCKING } from "../modules/reporting/summarize"

workflow SCREEN_WF {
    if (!params.receptor) { error "Parameter --receptor is required." }
    if (!params.library) { error "Parameter --library is required for screening." }

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

    ch_rec_pdb = PREPARE_RECEPTOR.out.prepared_pdb
    ch_ref_sdf = PREPARE_RECEPTOR.out.ref_native_sdf
    ch_rec_mae = PREPARE_RECEPTOR.out.prepared_mae
    ch_ref_mae = PREPARE_RECEPTOR.out.ref_native_mae

    // Prepare ligands dynamically (including native reference ligand)
    ch_lib = Channel.fromPath(file(params.library).toAbsolutePath(), checkIfExists: true)
        .map { lib -> [ [id: params.project_name, project: params.project_name], lib ] }

    PREPARE_LIGANDS(
        ch_lib,
        ch_ref_mae,
        params.dock_engine,
        params.batch_size,
        params.ligand_ph,
        params.ligand_ph_tol,
        params.ligprep_ff,
        params.max_states
    )

    // Parallel docking (Glide, GNINA, or Both)
        DOCK_LIBRARY(
            ch_rec_pdb,
            ch_ref_sdf,
            PREPARE_LIGANDS.out.chunks,
            ch_rec_mae,
            ch_ref_mae,
            PREPARE_LIGANDS.out.prepared_ligands,
            params.dock_engine,
            params.exhaustiveness,
            params.cnn_scoring,
            params.num_modes,
            params.autobox_add,
            params.top_k_frames,
            params.percentile_cutoff,
            params.membrane
        )
    // Summarize raw poses into consensus ranking CSV tables
    def bench_ligands_str = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""
    SUMMARIZE_DOCKING(
        params.project_name,
        DOCK_LIBRARY.out.summary_files,
        params.native_ligand,
        bench_ligands_str
    )
}
