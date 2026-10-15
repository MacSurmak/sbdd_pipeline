/*
 * subworkflows/prepare_receptor.nf
 * High-throughput parallel receptor preparation across Slurm cluster:
 * 1. EXTRACT_COMPLEX (Schrodinger): extracts pure protein PDB and full complex MAE (with cofactors & ligand).
 * 2. If membrane:
 *      IMMERS (Pure PPM 3.0 CLI): pure protein PDB -> oriented protein PDB + membrane_info.
 *      ORIENT_COMPLEX (Schrodinger): applies CA superposition matrix to the FULL complex structure,
 *          preserving all cofactors, bond orders, and relative native ligand geometry.
 *          Emits oriented receptor_raw.mae and ref_native.mae / ref_native.sdf.
 *    Else:
 *      Directly extracts receptor_raw.mae and ref_native.mae / ref_native.sdf.
 * 3. PREPWIZARD (Schrodinger): runs on receptor_raw.mae (protein + cofactors).
 * 4. If membrane:
 *      APPLY_MEMBRANE (Schrodinger Prime): injects implicit membrane slab properties for MM-GBSA.
 */

include { EXTRACT_COMPLEX } from "../modules/schrodinger/extract_complex"
include { IMMERS          } from "../modules/membrane/immers"
include { ORIENT_COMPLEX  } from "../modules/schrodinger/orient_complex"
include { PREPWIZARD      } from "../modules/schrodinger/prepwizard"
include { APPLY_MEMBRANE   } from "../modules/membrane/apply_membrane"
include { SPLIT_PREPARED   } from "../modules/schrodinger/split_prepared"

workflow PREPARE_RECEPTOR {
    take:
    ch_receptor       // channel: [ val(meta), path(complex_file) ]
    membrane_enabled  // val: true/false
    membrane_type     // val: "MOM"
    membrane_topo     // val: "in"
    membrane_chains   // val: "A"
    target_ph         // val: 7.4
    ph_tol            // val: 1.0
    min_rmsd          // val: 0.3
    cap_termini       // val: true

    main:
    EXTRACT_COMPLEX(
        ch_receptor,
        membrane_enabled
    )

    if (membrane_enabled.toString() == 'true') {
        IMMERS(
            EXTRACT_COMPLEX.out.pure_protein_pdb,
            membrane_type,
            membrane_topo,
            membrane_chains
        )

        ch_orient_in = EXTRACT_COMPLEX.out.pure_protein_pdb
            .join(IMMERS.out.oriented_pdb)
            .join(EXTRACT_COMPLEX.out.complex_mae)

        ORIENT_COMPLEX(
            ch_orient_in
        )

        ch_prep_in  = ORIENT_COMPLEX.out.oriented_mae
        ch_mem_json = IMMERS.out.membrane_json
        ch_mem_txt  = IMMERS.out.membrane_txt
    } else {
        ch_prep_in  = EXTRACT_COMPLEX.out.complex_mae
        ch_mem_json = Channel.empty()
        ch_mem_txt  = Channel.empty()
    }

    PREPWIZARD(
        ch_prep_in,
        target_ph,
        ph_tol,
        min_rmsd,
        cap_termini
    )

    if (membrane_enabled.toString() == 'true') {
        ch_apply_in = PREPWIZARD.out.prepared_mae.join(ch_mem_txt)
        APPLY_MEMBRANE(ch_apply_in)
        ch_final_complex_mae = APPLY_MEMBRANE.out.prepared_mae
    } else {
        ch_final_complex_mae = PREPWIZARD.out.prepared_mae
    }

    SPLIT_PREPARED(ch_final_complex_mae)

    // Emit structured 5-column receipt with absolute paths to canonical targets/<target>/frames/<frame>/ directory
    ch_rec_pair = SPLIT_PREPARED.out.prepared_mae.join(SPLIT_PREPARED.out.prepared_pdb)
    ch_lig_pair = SPLIT_PREPARED.out.ref_native_mae.join(SPLIT_PREPARED.out.ref_native_sdf)

    ch_receipt_data = ch_rec_pair
        .join(ch_lig_pair, remainder: true)
        .map { meta, r_mae, r_pdb, n_mae, n_sdf ->
            def frame_dir = file("${params.outdir}/targets/${meta.target}/frames/${meta.name}").toAbsolutePath()
            def sys_poses_dir = file("${params.outdir}/systems/${meta.system}/poses").toAbsolutePath()
            def mae_abs = "${frame_dir}/${r_mae.name}"
            def pdb_abs = "${frame_dir}/${r_pdb.name}"
            def n_mae_abs = (n_mae != null) ? "${sys_poses_dir}/${n_mae.name}" : "NONE"
            def n_sdf_abs = (n_sdf != null) ? "${sys_poses_dir}/${n_sdf.name}" : "NONE"
            "${meta.name}\t${mae_abs}\t${pdb_abs}\t${n_mae_abs}\t${n_sdf_abs}"
        }

    ch_receipt_data.collectFile(name: 'prep_receipt.tsv', storeDir: params.outdir, newLine: true)

    // Stream protonated in-situ poses directly to Antechamber for quantum parametrization into systems/<sys>/poses
    ch_poses_to_param = SPLIT_PREPARED.out.poses_manifest
        .flatMap { meta, json_file ->
            def records = new groovy.json.JsonSlurperClassic().parseText(json_file.text)
            records.findAll { !it.is_ion }.collect { rec ->
                def pose_meta = [
                    id: "${meta.name}_${rec.res_name}_pose",
                    name: "${meta.name}_${rec.res_name}_pose",
                    system: meta.system,
                    charge: rec.charge as int,
                    res_name: rec.res_name
                ]
                [ pose_meta, json_file.parent.resolve(rec.sdf) ]
            }
        }

    emit:
    prepared_mae     = SPLIT_PREPARED.out.prepared_mae
    ref_native_mae   = SPLIT_PREPARED.out.ref_native_mae
    ref_native_sdf   = SPLIT_PREPARED.out.ref_native_sdf
    membrane_json    = ch_mem_json
    poses_to_param   = ch_poses_to_param
}
