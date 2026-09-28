/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow with ensemble support
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'

workflow PREPARE_WF {
    if (!params.manifest) { error "FATAL: Parameter --manifest is required." }

    def manifestFile = file(params.manifest)
    if (!manifestFile.exists()) { error "FATAL: Manifest file not found at ${params.manifest}" }
    
    def yaml = new org.yaml.snakeyaml.Yaml()
    def manifestData = yaml.load(manifestFile.text)
    if (!manifestData.targets) { error "FATAL: 'targets' list missing in manifest.yaml" }

    ch_targets = Channel.fromList(manifestData.targets).map { target ->
        def target_id = target.id ?: "unknown"
        def struct_file = target.structure
        if (!struct_file) error "FATAL: Target ${target_id} has no 'structure' defined."
        
        def rec_path = file("${manifestFile.parent}/${struct_file}")
        if (!rec_path.exists()) {
            rec_path = file("${params.outdir}/${struct_file}")
        }
        if (!rec_path.exists()) error "FATAL: Structure file ${struct_file} not found for target ${target_id}."

        def meta = [
            id: target_id,
            project: params.project_name,
            native_ligand: target.native_ligand ?: "",
            cofactors: target.cofactors ?: [],
            is_prepared: false
        ]
        
        [ meta, rec_path.toAbsolutePath() ]
    }

    PREPARE_RECEPTOR(
        ch_targets,
        params.membrane,
        params.membrane_type,
        params.membrane_topo,
        params.membrane_chains,
        params.target_ph,
        params.ph_tolerance,
        params.minimize_rmsd,
        params.cap_termini
    )
}
