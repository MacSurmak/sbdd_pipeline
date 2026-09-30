/*
 * workflows/analyze.nf
 * Independent trajectory analysis workflow (PBC fixing, RMSD, MMPBSA)
 */

include { TRJCONV         } from '../modules/gromacs/trjconv'
include { ANALYZE_MD_TRAJ } from '../subworkflows/analyze_md_traj'
include { GMX_CLUSTER     } from '../modules/gromacs/cluster'

workflow ANALYZE_WF {

    // MODE 1: Clustering from DB plan
    if (params.cluster_plan) {
        def planFile = file(params.cluster_plan)
        def jsonSlurper = new groovy.json.JsonSlurper()
        def plan = jsonSlurper.parseText(planFile.text)

        def meta = [id: plan.traj_name, project: params.project_name]
        ch_in = Channel.value([ meta, file(plan.xtc), file(plan.tpr), file(plan.ref_pdb) ])
        
        GMX_CLUSTER(ch_in, plan.cutoff, plan.method)

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

    // MODE 2: Legacy manual analysis
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
        .combine(ch_topol)
        .combine(ch_ndx)
        .map { m, traj, tpr, top, ndx -> [m, traj, tpr, top, ndx] }

    TRJCONV(
        ch_trjconv_in,
        params.pbc_center,
        params.pbc_extract
    )

    ANALYZE_MD_TRAJ(
        TRJCONV.out.clean_traj,
        params.mmpbsa_start_ns,
        params.mmpbsa_step_frames,
        params.md_temperature,
        params.md_salt_conc,
        params.membrane,
        params.membrane_thickness,
        params.per_residue_decomp
    )
}
