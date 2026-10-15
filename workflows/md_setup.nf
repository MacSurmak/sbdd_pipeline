/*
 * workflows/md_setup.nf
 * Standalone setup routine to parametrize MD topologies
 */

include { SETUP_MD_SYSTEM } from '../subworkflows/setup_md_system'

process PUBLISH_MD_SYSTEM {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/systems/${meta.id}/md_topology" }, mode: 'copy'

    input:
    tuple val(meta), path(gro), path(top), path(ndx), path(posre)

    output:
    tuple val(meta), path("${meta.id}_step3.gro"), path("${meta.id}_topol.top"), path("${meta.id}_index.ndx"), path("posre_*.itp", includeInputs: true)

    script:
    """
    cp -L "${gro}" "${meta.id}_step3.gro"
    cp -L "${top}" "${meta.id}_topol.top"
    cp -L "${ndx}" "${meta.id}_index.ndx"
    """
}

workflow MD_SETUP_WF {
    if (!params.md_plan) { error "FATAL: Parameter --md_plan is required for md-setup." }

    def planFile = file(params.md_plan)
    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    ch_targets = Channel.fromList(planData).map { target ->
        def meta = [
            id: target.id,
            lig_resname: target.lig_resname ?: "LIG",
            extra_resnames: target.extra_resnames ?: [],
            project: params.project_name
        ]
        def rec_pdb = file(target.rec_pdb).toAbsolutePath()
        def lig_sdf = target.lig_sdf ? file(target.lig_sdf).toAbsolutePath() : []
        def extra_sdfs = target.extra_sdfs ? target.extra_sdfs.collect { file(it).toAbsolutePath() } : []
        def mol2s = target.mol2_paths.collect { file(it).toAbsolutePath() }
        def frcmods = target.frcmod_paths.collect { file(it).toAbsolutePath() }

        [ meta, rec_pdb, lig_sdf, extra_sdfs, mol2s, frcmods ]
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
        params.declash_tol,
        params.protein_forcefield,
        params.gaff_forcefield,
        params.lipid_forcefield
    )

    PUBLISH_MD_SYSTEM(SETUP_MD_SYSTEM.out.system_ready)

    PUBLISH_MD_SYSTEM.out.map { meta, gro, top, ndx, posre ->
        def sys_dir = file("${params.outdir}/systems/${meta.id}/md_topology").toAbsolutePath()
        "${meta.id}\t${sys_dir}/${gro.name}\t${sys_dir}/${top.name}\t${sys_dir}/${ndx.name}\n"
    }.collectFile(name: 'md_setup_receipt.tsv', storeDir: params.outdir)
}
