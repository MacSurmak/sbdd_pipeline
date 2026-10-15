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

        PREPARE_LIGANDS(
            ch_lig_in,
            params.batch_size ?: 15,
            params.ligand_ph ?: 7.4,
            params.ligand_ph_tol ?: 2.0,
            params.ligprep_ff ?: 16,
            params.max_states ?: 32
        )

        PREPARE_LIGANDS.out.prepared_ligands.map { meta, mae ->
            def base_id = meta.name ?: meta.id
            def stem = mae.name.replaceAll(/(\.maegz|\.mae)$/, "")
            "${base_id}\t${base_id}/${mae.name}\t${base_id}/${stem}.sdf"
        }.collectFile(name: 'prep_ligand_receipt.tsv', storeDir: params.outdir, newLine: true)
    }

    // 3. Process Poses for Quantum Parametrization
    def poseList = planData.poses ?: []
    if (poseList) {
        ch_poses_in = Channel.fromList(poseList).map { pose ->
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

        PARAMETRIZE_POSE(ch_poses_in)

        PARAMETRIZE_POSE.out.map { meta, mol2, frcmod ->
            "${meta.name}\tprepared_poses/${meta.name}/${mol2.name}\tprepared_poses/${meta.name}/${frcmod.name}"
        }.collectFile(name: 'prep_pose_receipt.tsv', storeDir: params.outdir, newLine: true)
    }
}
