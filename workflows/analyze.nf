/*
 * workflows/analyze.nf
 * Independent trajectory analysis workflow (PBC fixing, RMSD, MMPBSA)
 */

include { TRJCONV                   } from '../modules/gromacs/trjconv'
include { MD_METRICS                } from '../modules/gromacs/md_metrics'
include { GMX_MMPBSA                } from '../modules/gromacs/mmpbsa'
include { ANALYZE_MD_TRAJ           } from '../subworkflows/analyze_md_traj'
include { CLUSTER_ENSEMBLE          } from '../subworkflows/cluster_ensemble'
include { EXTRACT_MEDOID_COMPLEXES  } from '../modules/analysis/extract_medoids'

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
                ligand: item.ligand ?: "",
                clean_dir: item.clean_dir,
                project: params.project_name
            ]
            def xtc_files = item.xtcs ? item.xtcs.collect { file(it) } : [ file(item.xtc) ]
            def dry_top_f = item.dry_top ? file(item.dry_top) : []
            [ meta, xtc_files, file(item.tpr), dry_top_f ]
        }

        TRJCONV(
            ch_clean_in,
            params.pbc_center,
            "COMPLEX",
            params.water_names,
            params.ion_names,
            params.lipid_names
        )

        TRJCONV.out.clean_traj.map { meta, clean_xtc, ref_pdb, clean_ndx ->
            "${meta.id}\t${clean_xtc.name}\t${ref_pdb.name}\t${clean_ndx.name}"
        }.collectFile(name: 'clean_receipt.tsv', storeDir: params.outdir, newLine: true)

        return
    }

    // MODE 1.5: Distributed Slurm Extraction of Medoid Complexes
    if (params.extract_medoids_plan) {
        def planFile = file(params.extract_medoids_plan)
        def planList = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

        ch_medoids_in = Channel.fromList(planList).map { item ->
            def meta = [
                id: item.task_id,
                system: item.system_name,
                target: item.target_name,
                clustering: item.clustering_name,
                project: params.project_name
            ]
            def ref_pdb = file(item.ref_pdb)
            def xtcs = item.clean_xtcs.collect { file(it) }
            def spec_json = file("${workDir}/stage_${item.task_id}_spec.json")
            spec_json.text = groovy.json.JsonOutput.toJson(item)
            [ meta, ref_pdb, xtcs, spec_json ]
        }

        EXTRACT_MEDOID_COMPLEXES(ch_medoids_in)

        EXTRACT_MEDOID_COMPLEXES.out.receipt.collectFile(
            name: 'extract_medoids_receipt.tsv',
            storeDir: params.outdir,
            keepHeader: false
        )
        return
    }

    // MODE 2: Ensemble Clustering on Clean Trajectories (Analytical Subworkflow)
    if (params.cluster_plan) {
        def planFile = file(params.cluster_plan)
        def planList = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

        ch_cluster_in = Channel.fromList(planList).map { item ->
            def align_grp = item.align_group
            def cluster_grp = item.cluster_group
            def clean_cluster_tag = cluster_grp.replaceAll(/[^a-zA-Z0-9_]/, '_')
            def meta = [
                id: item.task_id,
                system: item.system,
                ligand: item.ligand ?: "",
                align_group: align_grp,
                cluster_group: cluster_grp,
                cluster_group_tag: clean_cluster_tag,
                replicate_names: item.replicate_names,
                project: params.project_name
            ]
            def xtcs = item.clean_xtcs.collect { file(it) }
            def ref_pdb = file(item.ref_pdb)
            def ndx_f = item.ndx ? file(item.ndx) : []
            [ meta, ref_pdb, xtcs, ndx_f ]
        }

        def sampling_start = params.analysis_start_ns ?: params.cluster_start_ns ?: 1.0
        def sampling_interval = params.analysis_sample_interval_ps ?: params.cluster_sample_interval_ps ?: 500.0

        CLUSTER_ENSEMBLE(
            ch_cluster_in,
            params.cluster_align_group ?: "Backbone",
            params.cluster_feature_group ?: "Ligand",
            params.cluster_dim_reduction ?: "tica",
            params.cluster_algorithm ?: "hdbscan",
            sampling_start,
            sampling_interval,
            params.cluster_tica_lag_ns ?: 2.0,
            params.cluster_variance_cutoff ?: 0.85,
            params.cluster_pocket_radius ?: 5.0,
            params.cluster_n_clusters ?: 3,
            params.cluster_rmsd_cutoff ?: 1.5,
            params.cluster_min_pop ?: 5.0,
            params.cluster_temp_k ?: 310.15
        )

        // Publish structured receipt mapping system and groups to published manifest
        CLUSTER_ENSEMBLE.out.manifest.map { meta, manifest_json ->
            def clean_tag = meta.cluster_group_tag ?: meta.cluster_group.replaceAll(/[^a-zA-Z0-9_]/, '_')
            def published_manifest = "${params.outdir}/reports/clusters/${meta.system}/${meta.align_group}_${clean_tag}/${manifest_json.name}"
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
            [ meta, file(item.clean_xtc), file(item.ref_pdb), file(item.dry_prmtop), file(item.clean_ndx) ]
        }

        ch_templates = file(params.mmpbsa_templates)

        def mmpbsa_start = params.analysis_start_ns ?: params.mmpbsa_start_ns ?: 1.0
        def mmpbsa_interval = params.analysis_sample_interval_ps ?: params.mmpbsa_sample_interval_ps ?: 500.0

        GMX_MMPBSA(
            ch_mmpbsa_in,
            ch_templates,
            mmpbsa_start,
            mmpbsa_interval,
            params.md_temperature,
            params.md_salt_conc,
            params.membrane,
            params.membrane_thickness,
            params.per_residue_decomp,
            params.mmpbsa_method ?: "auto",
            params.mmpbsa_gb_model ?: 2,
            params.mmpbsa_emem ?: 7.0,
            params.mmpbsa_indi ?: 4.0,
            params.mmpbsa_fillratio ?: 1.25
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
        .map { m, traj, tpr -> [m, traj, tpr, []] }

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
        file(params.mmpbsa_templates),
        params.mmpbsa_start_ns,
        params.mmpbsa_sample_interval_ps,
        params.md_temperature,
        params.md_salt_conc,
        params.membrane,
        params.membrane_thickness,
        params.per_residue_decomp,
        params.mmpbsa_method ?: "auto",
        params.mmpbsa_gb_model ?: 2,
        params.mmpbsa_emem ?: 7.0,
        params.mmpbsa_indi ?: 4.0,
        params.mmpbsa_fillratio ?: 1.25
    )
}
