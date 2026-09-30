/*
 * workflows/analyze.nf
 * Independent trajectory analysis workflow (PBC fixing, RMSD, MMPBSA)
 */

include { TRJCONV         } from '../modules/gromacs/trjconv'
include { MD_METRICS      } from '../modules/gromacs/md_metrics'
include { GMX_CLUSTER     } from '../modules/gromacs/cluster'
include { ANALYZE_MD_TRAJ } from '../subworkflows/analyze_md_traj'

workflow ANALYZE_WF {

    // MODE 1: Trajectory Cleaning (PBC whole -> center -> fit -> solvent strip)
    if (params.clean_plan) {
        def planFile = file(params.clean_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def planList = jsonSlurper.parseText(planFile.text)

        ch_clean_in = Channel.fromList(planList).map { item ->
            def meta = [id: item.name, rep: item.rep, project: params.project_name]
            [ meta, file(item.xtc), file(item.tpr) ]
        }

        TRJCONV(
            ch_clean_in,
            params.pbc_center,
            "COMPLEX"
        )

        // Emit receipt for database registration
        TRJCONV.out.clean_traj.map { meta, clean_xtc, ref_pdb ->
            "${meta.id}\t${clean_xtc.toAbsolutePath()}\t${ref_pdb.toAbsolutePath()}"
        }.collectFile(name: 'clean_receipt.tsv', storeDir: params.outdir, newLine: true)

        return
    }

    // MODE 2: Ensemble Clustering on Clean Trajectories
    if (params.cluster_plan) {
        def planFile = file(params.cluster_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def planList = jsonSlurper.parseText(planFile.text)

        ch_cluster_in = Channel.fromList(planList).map { item ->
            def meta = [id: item.task_id, project: params.project_name]
            def xtcs = item.clean_xtcs.collect { file(it) }
            def ref_pdb = file(item.ref_pdb)
            [ meta, xtcs, ref_pdb, item.cutoff, item.method ]
        }

        GMX_CLUSTER(ch_cluster_in)

        // Flatten output list of pdbs and emit receipt
        GMX_CLUSTER.out.cluster_pdbs.flatMap { m, pdbs ->
            def plist = pdbs instanceof List ? pdbs : [pdbs]
            plist.withIndex().collect { p, i ->
                def frame_id = "${m.id}_c${i+1}"
                "${frame_id}\t${p.toAbsolutePath()}"
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
            def meta = [id: item.name, rep: item.rep, project: params.project_name]
            [ meta, file(item.clean_xtc), file(item.ref_pdb) ]
        }

        MD_METRICS(ch_metrics_in)
        return
    }

    // MODE 4: Legacy manual analysis
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
