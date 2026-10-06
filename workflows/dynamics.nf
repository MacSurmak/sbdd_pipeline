/*
 * workflows/dynamics.nf
 * Standalone MD system setup, chunked simulation, and trajectory analysis
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { SETUP_MD_SYSTEM  } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD   } from '../subworkflows/run_chunked_md'
include { ANALYZE_MD_TRAJ  } from '../subworkflows/analyze_md_traj'

workflow DYNAMICS_WF {
    if (!params.md_protocol) { error "FATAL: Parameter --md_protocol is required. Did you run 'sbdd md-init'?" }

    if (params.md_plan) {
        def planFile = file(params.md_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def planData = jsonSlurper.parseText(planFile.text)

        ch_targets = Channel.fromList(planData).map { target ->
            def meta = [
                id: target.id,
                project: params.project_name,
                native_ligand: target.lig_name ?: "",
                replicates: target.reps_to_run ?: [1]
            ]
            def lig_f = target.lig_path ? file(target.lig_path).toAbsolutePath() : []
            [ meta, file(target.rec_path).toAbsolutePath(), lig_f, target.cofactors ?: [] ]
        }
    } else if (params.manifest) {
        def manifestFile = file(params.manifest)
        def yaml = new org.yaml.snakeyaml.Yaml()
        def manifestData = yaml.load(manifestFile.text)

        ch_targets = Channel.fromList(manifestData.targets).map { target ->
            def target_id = target.id
            def struct_file = target.structure
            def lig_file = target.ligand
            
            def rec_path = file("${params.outdir}/prepared/${target_id}_prepared.mae")
            if (!rec_path.exists()) rec_path = file("${params.outdir}/prepared/${target_id}_prepared.pdb")
            if (!rec_path.exists()) rec_path = file("${manifestFile.parent}/${struct_file}")
            if (!rec_path.exists()) rec_path = file("${params.outdir}/${struct_file}")
            
            def lig_path
            if (lig_file) {
                lig_path = file("${params.outdir}/${lig_file}")
            } else {
                lig_path = file("${params.outdir}/prepared/${target_id}_ref_native_prepared.sdf")
                if (!lig_path.exists()) lig_path = file("${params.outdir}/prepared/${target_id}_ref_native.sdf")
                if (!lig_path.exists()) lig_path = null
            }
            
            def meta = [
                id: target_id,
                project: params.project_name
            ]
            
            def lig_f = lig_path ? lig_path.toAbsolutePath() : []
            [ meta, rec_path.toAbsolutePath(), lig_f, target.cofactors ?: [] ]
        }
    } else {
        error "FATAL: Either --md_plan or --manifest must be provided for dynamics."
    }

    // 1. Omnivorous System Build
    SETUP_MD_SYSTEM(
        ch_targets,
        params.membrane,
        params.upper_lipids,
        params.upper_ratio,
        params.lower_lipids,
        params.lower_ratio,
        params.lipid_padding,
        params.box_padding,
        params.md_salt_conc,
        params.water_leaprc,
        params.water_box,
        params.box_type,
        params.ss_cutoff,
        params.water_names,
        params.ion_names,
        params.lipid_names,
        params.md_restraints_fc,
        params.declash_cutoff,
        params.declash_niter,
        params.declash_damping,
        params.declash_max_move,
        params.declash_jitter,
        params.declash_tol
    )

    ch_md_protocol = Channel.fromPath("${params.md_protocol}/*", checkIfExists: true).collect()

    // 2. Equilibration & Production using Protocol Playlist
    RUN_CHUNKED_MD(
        SETUP_MD_SYSTEM.out.system_ready,
        ch_md_protocol,
        params.md_total_ns,
        params.md_chunk_ns,
        params.md_replicates as int
    )
}
