/*
 * workflows/sbdd.nf
 * Complete end-to-end SBDD pipeline with multi-engine docking, rescoring,
 * Screening Funnel Triage (Physics + AI), Prime MM-GBSA, automated hit selection, and chunked GROMACS MD.
 */

include { STAGE_INPUTS      } from '../modules/common/stage_inputs'
include { PREPARE_RECEPTOR  } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS   } from '../subworkflows/prepare_ligands'
include { DOCK_LIBRARY      } from '../subworkflows/dock_library'
include { SETUP_MD_SYSTEM   } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD    } from '../subworkflows/run_chunked_md'
include { ANALYZE_MD_TRAJ   } from '../subworkflows/analyze_md_traj'
include { SUMMARIZE_DOCKING } from '../modules/reporting/summarize'
include { SELECT_HITS       } from '../modules/reporting/select_hits'
include { COUNT_MODELS      } from '../modules/common/count_models'

workflow SBDD_WF {
    if (!params.receptor) { error "Parameter --receptor is required." }

    def rec_path = file(params.receptor).toAbsolutePath().toString()
    def lib_path = params.library ? file(params.library).toAbsolutePath().toString() : 'NO_FILE'

    STAGE_INPUTS(rec_path, lib_path, 'NO_FILE')

    ch_rec_meta = STAGE_INPUTS.out.staged_receptor.map { rec ->
        [ [id: params.project_name, project: params.project_name], rec ]
    }

    // 1. Dynamic determination of NMR ensemble frames
    COUNT_MODELS(ch_rec_meta)

    ch_rec = STAGE_INPUTS.out.staged_receptor
        .combine(COUNT_MODELS.out.count_file)
        .flatMap { staged_rec, meta_cnt, count_txt ->
            def total = count_txt.text.trim() as int
            def req = params.model_index ? (params.model_index as int) : 0
            def models = (req > 0 && req <= total) ? [ req ] : (1..total).toList()
            models.collect { idx ->
                def frame_meta = [
                    id: "${params.project_name}_f${idx}",
                    frame: idx,
                    project: params.project_name
                ]
                [ frame_meta, staged_rec, idx ]
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

    // 2. Prepare Ligands & Docking
    if (params.library) {
        ch_lib = STAGE_INPUTS.out.staged_library
            .map { lib -> [ [id: params.project_name, project: params.project_name], lib ] }

        PREPARE_LIGANDS(
            ch_lib,
            PREPARE_RECEPTOR.out.ref_native_sdf,
            params.dock_engine ?: "both",
            params.batch_size ?: 15,
            params.ligand_ph ?: 7.4,
            params.ligand_ph_tol ?: 2.0,
            params.ligprep_ff ?: 16,
            params.max_states ?: 32
        )

        // 3. Multi-Engine Docking, Triage, and MM-GBSA
        DOCK_LIBRARY(
            PREPARE_RECEPTOR.out.prepared_pdb,
            PREPARE_RECEPTOR.out.ref_native_sdf,
            PREPARE_LIGANDS.out.chunks,
            PREPARE_RECEPTOR.out.prepared_mae,
            PREPARE_RECEPTOR.out.ref_native_mae,
            PREPARE_LIGANDS.out.prepared_ligands,
            params.dock_engine ?: "both",
            params.exhaustiveness ?: 16,
            params.cnn_scoring ?: "rescore",
            params.num_modes ?: 3,
            params.autobox_add ?: 4.0,
            params.top_k_frames ?: 3,
            params.percentile_cutoff ?: 0.50,
            params.membrane ?: false
        )

        // 4. Summarize results into master ranking table and cross-frame matrix
        def bench_ligands_sbdd = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""
        SUMMARIZE_DOCKING(
            params.project_name,
            DOCK_LIBRARY.out.summary_files,
            params.native_ligand ?: "",
            bench_ligands_sbdd
        )

        // 5. Automated Hit Selection for Molecular Dynamics
        ch_ranking_tuple = SUMMARIZE_DOCKING.out.ranking_csv
            .map { csv -> [ [id: params.project_name, project: params.project_name], csv ] }

        SELECT_HITS(
            ch_ranking_tuple,
            params.top_n_hits ?: 1,
            true
        )

        // 6. Downstream MD simulation
        if (params.run_md) {
            ch_hits = SELECT_HITS.out.hits_file
                .flatMap { meta, file -> file.readLines() }
                .map { line -> line.trim() }
                .filter { it.length() > 0 }
                .map { line ->
                    def parts = line.split(',')
                    def title = parts[0].trim()
                    def fid = parts.size() > 1 ? (parts[1].trim() as int) : 1
                    [ fid, title ]
                }

            ch_docked_by_frame = DOCK_LIBRARY.out.docked_poses
                .map { meta, pose -> [ (meta.frame as int), pose ] }

            ch_rec_by_frame = PREPARE_RECEPTOR.out.prepared_pdb
                .map { meta, pdb -> [ (meta.frame as int), pdb ] }

            ch_frame_data = ch_docked_by_frame.join(ch_rec_by_frame, by: 0)

            ch_md_targets = ch_hits.join(ch_frame_data, by: 0)
                .map { fid, hit_name, pose_file, rec_pdb ->
                    def target_id = hit_name.replaceAll("[^a-zA-Z0-9_-]", "_")
                    [ [id: target_id, frame: fid], pose_file, hit_name, rec_pdb ]
                }
            SETUP_MD_SYSTEM(
                ch_md_targets,
                params.membrane ?: false,
                params.upper_lipids ?: "POPC:POPE:CHL1",
                params.upper_ratio ?: "5:4:1",
                params.lower_lipids ?: "POPC:POPE:POPS:CHL1",
                params.lower_ratio ?: "4:4:1:1",
                params.lipid_padding ?: 10.0,
                params.box_padding ?: 15.0,
                params.md_salt_conc ?: 0.15,
                params.water_model ?: "opc"
            )

            RUN_CHUNKED_MD(
                SETUP_MD_SYSTEM.out.system_ready,
                params.membrane ?: false,
                params.md_temperature ?: 310.15,
                params.md_total_ns ?: 60,
                params.md_chunk_ns ?: 10,
                params.md_replicates ? (params.md_replicates as int) : 3,
                params.queue ?: "main",
                params.nvt_short_ps ?: 50,
                params.nvt_long_ps ?: 500,
                params.npt_berendsen_ps ?: 1000,
                params.npt_pr_ps ?: 1000
            )

            ANALYZE_MD_TRAJ(
                RUN_CHUNKED_MD.out.clean_trajectory,
                params.mmpbsa_cores ?: 16,
                params.mmpbsa_start_ns ?: 20.0,
                params.mmpbsa_step_frames ?: 40,
                params.md_temperature ?: 310.15,
                params.md_salt_conc ?: 0.15,
                params.membrane ?: false,
                params.membrane_thickness ?: 28.0,
                params.per_residue_decomp ?: true
            )
        }
    }
}
