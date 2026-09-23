/*
 * workflows/dynamics.nf
 * Standalone MD system setup, chunked simulation, and trajectory analysis
 */

include { STAGE_INPUTS     } from '../modules/common/stage_inputs'
include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { COUNT_MODELS     } from '../modules/common/count_models'
include { SETUP_MD_SYSTEM  } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD   } from '../subworkflows/run_chunked_md'
include { ANALYZE_MD_TRAJ  } from '../subworkflows/analyze_md_traj'

workflow DYNAMICS_WF {
    def rec_path = params.receptor ? file(params.receptor).toAbsolutePath().toString() : 'NO_FILE'
    def dock_path = params.docked_poses ? file(params.docked_poses).toAbsolutePath().toString() : 'NO_FILE'

    if (rec_path == 'NO_FILE') { error "No receptor found for dynamics. Provide --receptor." }

    STAGE_INPUTS(rec_path, 'NO_FILE', dock_path)

    def is_prepared = rec_path.endsWith('_prepared.pdb') || rec_path.endsWith('_prepared.mae')

    ch_rec_pdb = Channel.empty()
    ch_lig_sdf = Channel.empty()

    if (!is_prepared) {
        // Receptor is raw, we must PREPARE it first!
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
                    def frame_meta = [ id: "${params.project_name}_f${idx}", frame: idx, project: params.project_name ]
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
        ch_rec_pdb = PREPARE_RECEPTOR.out.prepared_pdb
        // If dock_path was not provided, use the native ligand we just extracted
        ch_lig_sdf = dock_path == 'NO_FILE' ? PREPARE_RECEPTOR.out.ref_native_sdf : STAGE_INPUTS.out.staged_docked.map { d -> [ [id: params.project_name, frame: 1], d ] }
    } else {
        // Already prepared, just use inputs directly
        ch_rec_pdb = STAGE_INPUTS.out.staged_receptor.map { r -> [ [id: params.project_name, frame: 1], r ] }
        ch_lig_sdf = STAGE_INPUTS.out.staged_docked.map { d -> [ [id: params.project_name, frame: 1], d ] }
    }

    // Parse target hit list
    def hits = params.md_hits ? params.md_hits.split(",").collect { it.trim() } : ["native"]

    ch_targets = Channel.fromList(hits)
        .combine(ch_lig_sdf)
        .combine(ch_rec_pdb, by: 0)
        .map { meta_frame, hit_name, lig_sdf, rec_pdb ->
            def target_id = hit_name.replaceAll("[^a-zA-Z0-9_-]", "_")
            [ [id: target_id, frame: meta_frame.frame], lig_sdf, hit_name, rec_pdb ]
        }

    // 1. Build systems with dynamic padding and water model
    SETUP_MD_SYSTEM(
        ch_targets,
        params.membrane.toString() == 'true',
        params.upper_lipids ?: "POPC:POPE:CHL1",            params.upper_ratio ?: "5:4:1",
            params.lower_lipids ?: "POPC:POPE:POPS:CHL1",
            params.lower_ratio ?: "4:4:1:1",
            params.lipid_padding ?: 10.0,
            params.box_padding ?: 15.0,
            params.md_salt_conc ?: 0.15,
            params.water_model ?: "opc"
        )

        // 2. Equilibration & dynamic atomic production MD
    RUN_CHUNKED_MD(
        SETUP_MD_SYSTEM.out.system_ready,
        params.membrane.toString() == 'true',
        params.md_temperature ?: 310.15,
        params.md_total_ns ?: 60,
        params.md_chunk_ns ?: 10,
        params.md_replicates ? (params.md_replicates as int) : 3,
        params.queue ?: "main",
        params.nvt_short_ps ?: 50,
        params.nvt_long_ps ?: 500,
        params.npt_berendsen_ps ?: 1000,
        params.npt_pr_ps ?: 1000
    )    // 3. Trajectory analysis & MM-PBSA
    ANALYZE_MD_TRAJ(
        RUN_CHUNKED_MD.out.clean_trajectory,
        params.mmpbsa_cores ?: 16,
        params.mmpbsa_start_ns ?: 20.0,
        params.mmpbsa_step_frames ?: 40,
        params.md_temperature ?: 310.15,
        params.md_salt_conc ?: 0.15,
        params.membrane.toString() == 'true',
        params.membrane_thickness ?: 28.0,
        params.per_residue_decomp.toString() == 'true'
    )
}
