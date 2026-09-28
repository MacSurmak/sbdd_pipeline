/*
 * workflows/dynamics.nf
 * Standalone MD system setup, chunked simulation, and trajectory analysis
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { SETUP_MD_SYSTEM  } from '../subworkflows/setup_md_system'
include { RUN_CHUNKED_MD   } from '../subworkflows/run_chunked_md'
include { ANALYZE_MD_TRAJ  } from '../subworkflows/analyze_md_traj'

workflow DYNAMICS_WF {
    if (!params.manifest) { error "FATAL: Parameter --manifest is required." }
    if (!params.md_protocol) { error "FATAL: Parameter --md_protocol is required. Did you run 'sbdd md-init'?" }

    def manifestFile = file(params.manifest)
    def yaml = new org.yaml.snakeyaml.Yaml()
    def manifestData = yaml.load(manifestFile.text)

    ch_targets = Channel.fromList(manifestData.targets).map { target ->
        def target_id = target.id
        def struct_file = target.structure
        def lig_file = target.ligand
        
        def rec_path = file("${launchDir}/prepared/${target_id}_prepared.pdb")
        if (!rec_path.exists()) rec_path = file("${manifestFile.parent}/${struct_file}")
        
        def lig_path
        if (lig_file) {
            lig_path = file("${launchDir}/${lig_file}")
        } else {
            lig_path = file("${launchDir}/prepared/${target_id}_ref_native.sdf")
            if (!lig_path.exists()) lig_path = file("NO_LIGAND")
        }
        
        def meta = [
            id: target_id,
            project: params.project_name
        ]
        
        [ meta, rec_path.toAbsolutePath(), lig_path.toAbsolutePath(), target.cofactors ?: [] ]
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
        params.water_model
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
