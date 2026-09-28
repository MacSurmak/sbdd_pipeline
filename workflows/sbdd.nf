/*
 * workflows/sbdd.nf
 * Complete end-to-end SBDD pipeline with multi-engine docking, rescoring,
 * Screening Funnel Triage (Physics + AI), Prime MM-GBSA, automated hit selection, and chunked GROMACS MD.
 */

include { PREPARE_RECEPTOR  } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS   } from '../subworkflows/prepare_ligands'
include { DOCK_LIBRARY      } from '../subworkflows/dock_library'
include { SETUP_MD_SYSTEM   } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD    } from '../subworkflows/run_chunked_md'
include { TRJCONV           } from '../modules/gromacs/trjconv'
include { ANALYZE_MD_TRAJ   } from '../subworkflows/analyze_md_traj'
include { SUMMARIZE_DOCKING } from '../modules/reporting/summarize'
include { SELECT_HITS       } from '../modules/reporting/select_hits'

workflow SBDD_WF {
    if (!params.receptor) { error "Parameter --receptor is required." }

    ch_rec = Channel.fromPath(file(params.receptor).toAbsolutePath(), checkIfExists: true)
        .map { rec -> [ [id: params.project_name, project: params.project_name], rec ] }

    // 1. Receptor Preparation
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

    // 2. Prepare Ligands & Docking
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

        // 3. Multi-Engine Docking, Triage, and MM-GBSA
        DOCK_LIBRARY(
            PREPARE_RECEPTOR.out.prepared_pdb,
            PREPARE_RECEPTOR.out.ref_native_sdf,
            PREPARE_LIGANDS.out.chunks,
            PREPARE_RECEPTOR.out.prepared_mae,
            PREPARE_RECEPTOR.out.ref_native_mae,
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

        // 4. Summarize results into master ranking table and cross-frame matrix
        def bench_ligands_sbdd = params.benchmark_ligands ? (params.benchmark_ligands instanceof List ? params.benchmark_ligands.join(" ") : params.benchmark_ligands) : ""
        SUMMARIZE_DOCKING(
            params.project_name,
            DOCK_LIBRARY.out.summary_files,
            params.native_ligand,
            bench_ligands_sbdd
        )

        // 5. Automated Hit Selection for Molecular Dynamics
        ch_ranking_tuple = SUMMARIZE_DOCKING.out.ranking_csv
            .map { csv -> [ [id: params.project_name, project: params.project_name], csv ] }

        SELECT_HITS(
            ch_ranking_tuple,
            params.top_n_hits,
            true,
            params.native_ligand
        )

        // 6. Downstream MD simulation
        if (params.run_md.toString() == 'true') {
            if (!params.nvt_short_ps || !params.nvt_long_ps || !params.npt_berendsen_ps || !params.npt_pr_ps) {
                error "FATAL: MD equilibration parameters (--nvt_short_ps, --nvt_long_ps, --npt_berendsen_ps, --npt_pr_ps) are not set! Check md.yaml."
            }
            if (!params.md_total_ns || !params.md_chunk_ns || !params.md_temperature || !params.md_replicates) {
                error "FATAL: MD production parameters (--md_total_ns, --md_chunk_ns, --md_temperature, --md_replicates) are not set! Check md.yaml."
            }
            if (!params.water_model || !params.box_padding || !params.md_salt_conc) {
                error "FATAL: MD system build parameters (--water_model, --box_padding, --md_salt_conc) are not set! Check md.yaml."
            }

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
                params.membrane,
                params.upper_lipids,
                params.upper_ratio,
                params.lower_lipids,
                params.lower_ratio,
                params.lipid_padding,
                params.box_padding,
                params.md_salt_conc,
                params.water_model
            )

        ch_mdp_templates = Channel.fromPath("${projectDir}/assets/mdp/*.mdp", checkIfExists: true).collect()

        RUN_CHUNKED_MD(
            SETUP_MD_SYSTEM.out.system_ready,
            ch_mdp_templates,
            params.membrane,
            params.md_temperature,
            params.md_total_ns,
            params.md_chunk_ns,
            params.md_replicates as int,
            params.queue,
            params.nvt_short_ps,
            params.nvt_long_ps,
            params.npt_berendsen_ps,
            params.npt_pr_ps
        )

        // 7. Analysis (PBC fixing, RMSD, MMPBSA) strictly on completed raw trajectories
        ch_trjconv_in = RUN_CHUNKED_MD.out.raw_trajectory
            .join(RUN_CHUNKED_MD.out.sys_files, by: 0)

        TRJCONV(
            ch_trjconv_in,
            params.pbc_center,
            params.pbc_extract
        )

        ANALYZE_MD_TRAJ(
            TRJCONV.out.clean_traj,
            params.mmpbsa_start_ns,
            params.mmpbsa_step_frames,
            params.md_temperature,
            params.md_salt_conc,
            params.membrane,
            params.membrane_thickness,
            params.per_residue_decomp
        )
        }
    }
}
