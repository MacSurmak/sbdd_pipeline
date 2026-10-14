/*
 * workflows/md_setup.nf
 * Standalone setup routine to parametrize MD topologies
 */

include { SETUP_MD_SYSTEM } from '../subworkflows/setup_md_system'

process PUBLISH_MD_SYSTEM {
    tag "${meta.id}"
    publishDir "${params.outdir}/md_systems/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(gro), path(top), path(ndx), path(posre)

    output:
    tuple val(meta), path(gro), path(top), path(ndx), path(posre)

    script:
    """
    # Publish passthrough
    """
}

workflow MD_SETUP_WF {
    if (!params.md_plan) { error "FATAL: Parameter --md_plan is required for md-setup." }

    def planFile = file(params.md_plan)
    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    ch_targets = Channel.fromList(planData).map { target ->
        def meta = [
            id: target.id,
            project: params.project_name,
            native_ligand: target.lig_name ?: ""
        ]
        def lig_f = target.lig_path ? file(target.lig_path).toAbsolutePath() : []
        [ meta, file(target.rec_path).toAbsolutePath(), lig_f, target.cofactors ?: [] ]
    }

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

    PUBLISH_MD_SYSTEM(SETUP_MD_SYSTEM.out.system_ready)

    PUBLISH_MD_SYSTEM.out.map { meta, gro, top, ndx, posre ->
        "${meta.id}\t${params.outdir}/md_systems/${meta.id}/${gro.name}\t${params.outdir}/md_systems/${meta.id}/${top.name}\t${params.outdir}/md_systems/${meta.id}/${ndx.name}\n"
    }.collectFile(name: 'md_setup_receipt.tsv', storeDir: params.outdir)
}
