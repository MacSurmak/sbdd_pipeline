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
        params.model_index ?: 0,
        params.native_ligand,
        params.keep_cofactors ?: "",
        params.membrane.toString() == 'true',
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

    // Prepare ligands dynamically (including native reference ligand)
    ch_lib = Channel.fromPath(file(params.library).toAbsolutePath(), checkIfExists: true)
        .map { lib -> [ [id: params.project_name, project: params.project_name], lib ] }

    PREPARE_LIGANDS(
        ch_lib,
        ch_ref_mae,
        params.dock_engine ?: "both",
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
            params.dock_engine ?: "both",
            params.exhaustiveness ?: 16,
            params.cnn_scoring ?: "rescore",
            params.num_modes ?: 3,
            params.autobox_add ?: 4.0,
            params.top_k_frames ?: 3,
            params.percentile_cutoff ?: 0.50,
            params.membrane.toString() == 'true'
        )
    // Summarize raw poses into consensus ranking CSV tables
    def bench_ligands_str = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""
    SUMMARIZE_DOCKING(
        params.project_name,
        DOCK_LIBRARY.out.summary_files,
        params.native_ligand ?: "",
        bench_ligands_str
    )
}
