/*
 * workflows/analyze.nf
 * Independent trajectory analysis workflow (PBC fixing, RMSD, MMPBSA)
 */

include { TRJCONV          } from '../modules/gromacs/trjconv'
include { MD_METRICS       } from '../modules/gromacs/md_metrics'
include { GMX_MMPBSA       } from '../modules/gromacs/mmpbsa'
include { ANALYZE_MD_TRAJ  } from '../subworkflows/analyze_md_traj'
include { CLUSTER_ENSEMBLE } from '../subworkflows/cluster_ensemble'

workflow ANALYZE_WF {

    // MODE 1: Trajectory Cleaning (PBC whole -> center -> fit -> solvent strip)
    if (params.clean_plan) {
        def planFile = file(params.clean_plan)
        def planList = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

        ch_clean_in = Channel.fromList(planList).map { item ->
            def meta = [
                id: item.name,
                rec_id: item.rec_id,
                rep: item.rep,
                clean_dir: item.clean_dir,
                project: params.project_name
            ]
            def xtc_files = item.xtcs ? item.xtcs.collect { file(it) } : [ file(item.xtc) ]
            [ meta, xtc_files, file(item.tpr) ]
        }

        TRJCONV(
            ch_clean_in,
            params.pbc_center,
            "COMPLEX",
            params.water_names,
            params.ion_names,
            params.lipid_names
        )

        TRJCONV.out.clean_traj.map { meta, clean_xtc, ref_pdb ->
            "${meta.id}\t${clean_xtc.name}\t${ref_pdb.name}"
        }.collectFile(name: 'clean_receipt.tsv', storeDir: params.outdir, newLine: true)

        return
    }

    // MODE 2: Ensemble Clustering on Clean Trajectories (4-Stage Atomic Subworkflow)
    if (params.cluster_plan) {
        def planFile = file(params.cluster_plan)
        def planList = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

        ch_cluster_in = Channel.fromList(planList).map { item ->
            def align_grp = item.align_group
            def cluster_grp = item.cluster_group
            def meta = [
                id: item.task_id,
                system: item.system,
                ligand: item.ligand ?: "",
                align_group: align_grp,
                cluster_group: cluster_grp,
                replicate_names: item.replicate_names,
                project: params.project_name
            ]
            def xtcs = item.clean_xtcs.collect { file(it) }
            def ref_pdb = file(item.ref_pdb)
            def ndx_f = item.ndx ? file(item.ndx) : file("${workDir}/NO_FILE")
            def tmpl_sdf = item.template_sdf ? file(item.template_sdf) : file("${workDir}/NO_FILE")
            [ meta, ref_pdb, xtcs, ndx_f, tmpl_sdf ]
        }

        CLUSTER_ENSEMBLE(
            ch_cluster_in,
            params.cluster_align_group ?: "Backbone",
            params.cluster_feature_group ?: "Pocket",
            params.cluster_dim_reduction ?: "pca",
            params.cluster_tica_lag_ns ?: 2.0,
            params.cluster_variance_cutoff ?: 0.85,
            params.cluster_pocket_radius ?: 5.0,
            params.cluster_n_clusters ?: 3,
            params.cluster_temp_k ?: 310.15
        )

        // Publish structured receipt mapping system and groups to published manifest
        CLUSTER_ENSEMBLE.out.manifest.map { meta, manifest_json ->
            def published_manifest = "${params.outdir}/reports/clusters/${meta.system}/${meta.align_group}_${meta.cluster_group}/${manifest_json.name}"
            "${meta.system}\t${meta.align_group}\t${meta.cluster_group}\t${published_manifest}"
        }.collectFile(name: 'cluster_receipt.tsv', storeDir: params.outdir, newLine: true)

        return
    }

    // MODE 3: Trajectory Metrics (RMSD, RMSF, Contacts)
    if (params.metrics_plan) {
        def planFile = file(params.metrics_plan)
        def planList = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

        ch_metrics_in = Channel.fromList(planList).map { item ->
            def meta = [
                id: item.name,
                rec_id: item.rec_id,
                rep: item.rep,
                ligand: item.ligand ?: "",
                cofactors: item.cofactors ?: [],
                metrics_dir: item.metrics_dir,
                project: params.project_name
            ]
            [ meta, file(item.clean_xtc), file(item.ref_pdb) ]
        }

        MD_METRICS(ch_metrics_in, params.metric_heavy, params.metric_hyd, params.metric_pol)
        return
    }

    // MODE 4: Standalone MM-PBSA from DB
    if (params.mmpbsa_plan) {
        def planFile = file(params.mmpbsa_plan)
        def planList = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

        ch_mmpbsa_in = Channel.fromList(planList).map { item ->
            def meta = [
                id: item.name,
                rec_id: item.rec_id,
                rep: item.rep,
                ligand: item.ligand ?: "",
                project: params.project_name
            ]
            [ meta, file(item.clean_xtc), file(item.ref_pdb), file(item.topol), file(item.ndx) ]
        }

        GMX_MMPBSA(
            ch_mmpbsa_in,
            params.mmpbsa_start_ns,
            params.mmpbsa_step_frames,
            params.md_temperature,
            params.md_salt_conc,
            params.membrane,
            params.membrane_thickness,
            params.per_residue_decomp
        )
        return
    }

    // MODE 5: Legacy manual analysis
    if (!params.traj || !params.tpr || !params.topol || !params.ndx) {
        error "Analysis workflow requires explicitly provided --traj, --tpr, --topol, and --ndx parameters."
    }

    ch_traj  = Channel.fromPath(params.traj, checkIfExists: true)
    ch_tpr   = Channel.fromPath(params.tpr, checkIfExists: true)
    ch_topol = Channel.fromPath(params.topol, checkIfExists: true)
    ch_ndx   = Channel.fromPath(params.ndx, checkIfExists: true)

    // Dummy meta block since we operate on single trajectories in standalone mode
    def meta = [id: params.project_name ?: "md_analysis", rep: 1, ligand: params.native_ligand ?: ""]

    ch_trjconv_in = Channel.value(meta)
        .combine(ch_traj)
        .combine(ch_tpr)
        .map { m, traj, tpr -> [m, traj, tpr] }

    TRJCONV(
        ch_trjconv_in,
        params.pbc_center,
        params.pbc_extract,
        params.water_names,
        params.ion_names,
        params.lipid_names
    )

    ch_analyze_in = TRJCONV.out.clean_traj
        .combine(ch_topol)
        .combine(ch_ndx)
        .map { m, clean_xtc, ref_pdb, top, ndx -> [m, clean_xtc, ref_pdb, top, ndx] }

    ANALYZE_MD_TRAJ(
        ch_analyze_in,
        params.mmpbsa_start_ns,
        params.mmpbsa_step_frames,
        params.md_temperature,
        params.md_salt_conc,
        params.membrane,
        params.membrane_thickness,
        params.per_residue_decomp
    )
}
