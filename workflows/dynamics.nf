/*
 * workflows/dynamics.nf
 * Standalone MD system setup, chunked simulation, and trajectory analysis
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { SETUP_MD_SYSTEM  } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD   } from '../subworkflows/run_chunked_md'
include { ANALYZE_MD_TRAJ  } from '../subworkflows/analyze_md_traj'

workflow DYNAMICS_WF {
    if (!params.receptor) { error "Parameter --receptor is required for dynamics." }

    ch_rec = Channel.fromPath(file(params.receptor).toAbsolutePath(), checkIfExists: true)
        .map { rec -> [ [id: params.project_name, project: params.project_name], rec ] }

    // We completely trust native Nextflow caching (-resume). Always pass through PREPARE_RECEPTOR.
    PREPARE_RECEPTOR(
        ch_rec,
        params.model_index ?: 0,
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

    if (params.docked_poses) {
        ch_lig_sdf = Channel.fromPath(file(params.docked_poses).toAbsolutePath(), checkIfExists: true)
            .map { d -> [ [id: params.project_name, frame: 1], d ] }
    } else {
        ch_lig_sdf = PREPARE_RECEPTOR.out.ref_native_sdf
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

        ch_mdp_templates = Channel.fromPath("${projectDir}/assets/mdp/*.mdp", checkIfExists: true).collect()

    // 2. Equilibration & dynamic atomic production MD
    if (!params.nvt_short_ps || !params.nvt_long_ps || !params.npt_berendsen_ps || !params.npt_pr_ps) {
        error "FATAL: MD equilibration parameters (--nvt_short_ps, --nvt_long_ps, --npt_berendsen_ps, --npt_pr_ps) are not set! Check md.yaml."
    }
    if (!params.md_total_ns || !params.md_chunk_ns || !params.md_temperature || !params.md_replicates) {
        error "FATAL: MD production parameters (--md_total_ns, --md_chunk_ns, --md_temperature, --md_replicates) are not set! Check md.yaml."
    }

    RUN_CHUNKED_MD(
        SETUP_MD_SYSTEM.out.system_ready,
        ch_mdp_templates,
        params.membrane.toString() == 'true',
        params.md_temperature,
        params.md_total_ns,
        params.md_chunk_ns,
        params.md_replicates as int,
        params.queue,
        params.nvt_short_ps,
        params.nvt_long_ps,
        params.npt_berendsen_ps,
        params.npt_pr_ps
    )    // Note: Analysis is intentionally decoupled from MD generation.
    // Use `sbdd md-analyze` to process the RAW concatenated trajectories.
}
