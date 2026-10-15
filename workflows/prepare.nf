/*
 * workflows/prepare.nf
 * Standalone receptor and ligand preparation workflow with ensemble support
 */

include { PREPARE_RECEPTOR } from '../subworkflows/prepare_receptor'
include { PREPARE_LIGANDS  } from '../subworkflows/prepare_ligands'
include { PARAMETRIZE_POSE } from '../modules/amber/antechamber'

workflow PREPARE_WF {
    if (!params.plan) { error "FATAL: Parameter --plan (JSON) is required for database-driven preparation." }
    
    def planFile = file(params.plan)
    if (!planFile.exists()) { error "FATAL: Plan file not found at ${params.plan}" }
    
    def planData = new groovy.json.JsonSlurperClassic().parseText(planFile.text)

    // 1. Process Receptors if present in plan
    def receptorList = planData.receptors ?: (planData instanceof List ? planData : [])
    ch_extracted_poses = Channel.empty()

    if (receptorList) {
        ch_targets = Channel.fromList(receptorList).map { target ->
            def rec_path = file(target.path)
            if (!rec_path.exists()) {
                error "FATAL: Path explicitly provided by DB does not exist: ${target.path}"
            }

            def meta = [
                id: target.name,
                name: target.name,
                prep_mode: target.prep_mode,
                project: params.project_name,
                native_ligand: target.native_ligand ?: "",
                cofactors: target.cofactors ?: []
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

        ch_extracted_poses = PREPARE_RECEPTOR.out.poses_to_param
    }

    // 2. Process Ligands if present in plan
    def ligandList = planData.ligands ?: []
    if (ligandList) {
        ch_lig_in = Channel.fromList(ligandList).map { lig ->
            def lig_path = file(lig.path)
            if (!lig_path.exists()) {
                error "FATAL: Ligand file path does not exist: ${lig.path}"
            }
            def meta = [
                id: lig.name,
                name: lig.name,
                project: params.project_name,
                status: lig.status ?: 'unprepared',
                is_native: lig.is_native ?: false
            ]
            [ meta, lig_path.toAbsolutePath() ]
        }

        if (!params.batch_size || !params.ligand_ph || !params.ligand_ph_tol || !params.ligprep_ff || !params.max_states) {
            error "FATAL: Ligand preparation parameters missing in configuration! Check project.yaml."
        }

        PREPARE_LIGANDS(
            ch_lig_in,
            params.batch_size,
            params.ligand_ph,
            params.ligand_ph_tol,
            params.ligprep_ff,
            params.max_states
        )

        PREPARE_LIGANDS.out.prepared_ligands.map { meta, mae ->
            def base_id = meta.name ?: meta.id
            def stem = mae.name.replaceAll(/(\.maegz|\.mae)$/, "")
            def lig_dir = file("${params.outdir}/ligands/${base_id}").toAbsolutePath()
            "${base_id}\t${lig_dir}/${mae.name}\t${lig_dir}/${stem}.sdf"
        }.collectFile(name: 'prep_ligand_receipt.tsv', storeDir: params.outdir, newLine: true)
    }

    // 3. Process Poses for Quantum Parametrization (Standalone Poses + In-Situ Poses from PrepWizard)
    def poseList = planData.poses ?: []
    ch_plan_poses = Channel.empty()
    if (poseList) {
        ch_plan_poses = Channel.fromList(poseList).map { pose ->
            def sdf_path = file(pose.path)
            if (!sdf_path.exists()) {
                error "FATAL: Pose SDF path explicitly provided by DB does not exist: ${pose.path}"
            }
            def meta = [
                id: pose.name,
                name: pose.name,
                charge: pose.charge,
                res_name: pose.res_name ?: "LIG"
            ]
            [ meta, sdf_path.toAbsolutePath() ]
        }
    }

    ch_all_poses = ch_plan_poses.mix(ch_extracted_poses)

    PARAMETRIZE_POSE(ch_all_poses)

    PARAMETRIZE_POSE.out.map { meta, mol2, frcmod ->
        def pose_dir = file("${params.outdir}/prepared_poses/${meta.name}").toAbsolutePath()
        "${meta.name}\t${pose_dir}/${mol2.name}\t${pose_dir}/${frcmod.name}"
    }.collectFile(name: 'prep_pose_receipt.tsv', storeDir: params.outdir, newLine: true)
}
