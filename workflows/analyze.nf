/*
 * workflows/analyze.nf
 * Independent trajectory analysis workflow (PBC fixing, RMSD, MMPBSA)
 */

include { TRJCONV         } from '../modules/gromacs/trjconv'
include { MD_METRICS      } from '../modules/gromacs/md_metrics'
include { GMX_CLUSTER     } from '../modules/gromacs/cluster'
include { GMX_MMPBSA      } from '../modules/gromacs/mmpbsa'
include { ANALYZE_MD_TRAJ } from '../subworkflows/analyze_md_traj'

workflow ANALYZE_WF {

    // MODE 1: Trajectory Cleaning (PBC whole -> center -> fit -> solvent strip)
    if (params.clean_plan) {
        def planFile = file(params.clean_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def planList = jsonSlurper.parseText(planFile.text)

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
            "COMPLEX"
        )

        TRJCONV.out.clean_traj.map { meta, clean_xtc, ref_pdb ->
            "${meta.id}\t${clean_xtc.name}\t${ref_pdb.name}"
        }.collectFile(name: 'clean_receipt.tsv', storeDir: params.outdir, newLine: true)

        return
    }

    // MODE 2: Ensemble Clustering on Clean Trajectories
    if (params.cluster_plan) {
        def planFile = file(params.cluster_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def planList = jsonSlurper.parseText(planFile.text)

        ch_cluster_in = Channel.fromList(planList).map { item ->
            def meta = [
                id: item.task_id,
                ligand: item.ligand ?: "",
                max_clusters: item.max_clusters ?: 5,
                project: params.project_name
            ]
            def xtcs = item.clean_xtcs.collect { file(it) }
            def ref_pdb = file(item.ref_pdb)
            [ meta, xtcs, ref_pdb, item.cutoff, item.method, item.dt, item.group ]
        }

        GMX_CLUSTER(ch_cluster_in)

        GMX_CLUSTER.out.cluster_pdbs.flatMap { m, pdbs ->
            def plist = pdbs instanceof List ? pdbs : [pdbs]
            plist.sort { it.name }.collect { p ->
                def frame_id = p.baseName
                def idx_str = frame_id.tokenize('_').last()
                def pose_file = file("${p.parent}/${m.id}_pose_${idx_str}.pdb")
                def pose_path = pose_file.exists() ? pose_file.toAbsolutePath() : "NO_POSE"
                "${frame_id}\t${p.toAbsolutePath()}\t${pose_path}"
            }
        }.collectFile(name: 'cluster_receipt.tsv', storeDir: params.outdir, newLine: true)

        return
    }

    // MODE 3: Trajectory Metrics (RMSD, RMSF, Contacts)
    if (params.metrics_plan) {
        def planFile = file(params.metrics_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def planList = jsonSlurper.parseText(planFile.text)

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
            [ meta, file(item.clean_xtc), file(item.ref_gro) ]
        }

        MD_METRICS(ch_metrics_in)
        return
    }

    // MODE 4: Standalone MM-PBSA from DB
    if (params.mmpbsa_plan) {
        def planFile = file(params.mmpbsa_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def planList = jsonSlurper.parseText(planFile.text)

        ch_mmpbsa_in = Channel.fromList(planList).map { item ->
            def meta = [
                id: item.name,
                rec_id: item.rec_id,
                rep: item.rep,
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
    def meta = [id: params.project_name ?: "md_analysis", rep: 1]

    ch_trjconv_in = Channel.value(meta)
        .combine(ch_traj)
        .combine(ch_tpr)
        .map { m, traj, tpr -> [m, traj, tpr] }

    TRJCONV(
        ch_trjconv_in,
        params.pbc_center,
        params.pbc_extract
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
