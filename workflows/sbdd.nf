/*
 * workflows/sbdd.nf
 * Complete end-to-end SBDD pipeline
 * Complete end-to-end SBDD pipeline with multi-engine docking, rescoring,
 * Prime MM-GBSA, automated hit selection, and chunked GROMACS MD.
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

workflow SBDD_WF {
    def rec_file = file(params.receptor)
    def lib_file = params.library ? file(params.library) : file('NO_FILE')

    STAGE_INPUTS(rec_file, lib_file, file('NO_FILE'))

    // 1. Prepare Receptor
    ch_rec = STAGE_INPUTS.out.staged_receptor
        .map { r -> [ [id: params.project_name], r, params.model_index ?: 1 ] }

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

    // 2. Prepare Ligands & Docking
    if (params.library) {
        ch_lib = STAGE_INPUTS.out.staged_library
            .map { l -> [ [id: params.project_name], l ] }

        PREPARE_LIGANDS(
            ch_lib,
            params.dock_engine ?: "gnina",
            params.batch_size ?: 15,
            7.4,
            2.0,
            16,
            32
        )

        // 3. Multi-Engine Docking & MM-GBSA
        DOCK_LIBRARY(
            PREPARE_RECEPTOR.out.prepared_pdb,
            PREPARE_RECEPTOR.out.ref_native_sdf,
            PREPARE_LIGANDS.out.chunks,
            PREPARE_RECEPTOR.out.prepared_mae,
            PREPARE_RECEPTOR.out.ref_native_mae,
            PREPARE_LIGANDS.out.prepared_ligands,
            params.dock_engine ?: "gnina",
            params.exhaustiveness ?: 16,
            params.cnn_scoring ?: "rescore",
            params.num_modes ?: 3,
            params.autobox_add ?: 4.0
        )

        // 4. Summarize results into master ranking table
        SUMMARIZE_DOCKING(
            DOCK_LIBRARY.out.summary_bundle
        )

        // 5. Automated Hit Selection for Molecular Dynamics
        SELECT_HITS(
            SUMMARIZE_DOCKING.out.ranking_csv,
            params.top_n_hits ?: 1,
            true
        )

        // 6. Downstream MD simulation
        if (params.run_md) {
            def mdp_minim = file("${projectDir}/assets/mdp/minim.mdp")
            def mdp_nvt   = file("${projectDir}/assets/mdp/nvt.mdp")
            def mdp_npt   = file("${projectDir}/assets/mdp/npt.mdp")
            def mdp_md    = file("${projectDir}/assets/mdp/md.mdp")

            ch_md_targets = SELECT_HITS.out.hits_file
                .splitText()
                .flatMap { meta, file -> file.readLines() }
                .map { line -> line.trim() }
                .filter { it.length() > 0 }
                .combine(DOCK_LIBRARY.out.docked_poses)
                .combine(PREPARE_RECEPTOR.out.prepared_pdb)
                .map { hit_name, meta_dock, pose_file, meta_rec, rec_pdb ->
                    def target_id = hit_name.replaceAll("[^a-zA-Z0-9_-]", "_")
                    [ [id: target_id], pose_file, hit_name, rec_pdb ]
                }

            SETUP_MD_SYSTEM(
                ch_md_targets,
                params.membrane ?: false,
                params.upper_lipids ?: "POPC:POPE:CHL1",
                params.upper_ratio ?: "5:4:1",
                params.lower_lipids ?: "POPC:POPE:POPS:CHL1",
                params.lower_ratio ?: "4:4:1:1",
                10.0,
                15.0,
                params.md_salt_conc ?: 0.15,
                "opc"
            )

            RUN_CHUNKED_MD(
                SETUP_MD_SYSTEM.out.system_ready,
                mdp_minim,
                mdp_nvt,
                mdp_npt,
                mdp_md,
                params.membrane ?: false,
                params.md_temperature ?: 310.15,
                params.md_total_ns ?: 60,
                params.md_chunk_ns ?: 10,
                params.md_replicates ?: 3
            )

            ANALYZE_MD_TRAJ(
                RUN_CHUNKED_MD.out.clean_trajectory,
                params.mmpbsa_cores ?: 16,
                params.mmpbsa_start_ns ?: 20.0,
                params.mmpbsa_step_frames ?: 40,
                params.md_temperature ?: 310.15,
                true
            )
        }
    }
}
