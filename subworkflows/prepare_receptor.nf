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

    // Emit structured 3-column receipt strictly for ReceptorFrame
    ch_receipt_data = SPLIT_PREPARED.out.prepared_mae
        .join(SPLIT_PREPARED.out.prepared_pdb)
        .map { meta, r_mae, r_pdb ->
            def frame_dir = file("${params.outdir}/targets/${meta.target}/frames/${meta.name}").toAbsolutePath()
            def mae_abs = "${frame_dir}/${r_mae.name}"
            def pdb_abs = "${frame_dir}/${r_pdb.name}"
            "${meta.name}\t${mae_abs}\t${pdb_abs}"
        }

    ch_receipt_data.collectFile(name: 'prep_receipt.tsv', storeDir: params.outdir, newLine: true)

    // Stream protonated in-situ poses directly to Antechamber for quantum parametrization into substances/<inchikey>/poses/<pose_name>
    ch_poses_to_param = SPLIT_PREPARED.out.poses_manifest
        .flatMap { meta, json_file ->
            def records = new groovy.json.JsonSlurperClassic().parseText(json_file.text)
            def het_map = meta.heterogens ? meta.heterogens.collectEntries { [it.res_name, it] } : [:]

            records.collect { rec ->
                def het_info = het_map[rec.res_name]
                def p_name = het_info ? het_info.pose_name : "${meta.name}_${rec.res_name}_crystal"
                def p_inchikey = het_info ? het_info.inchikey : "UNKNOWN"
                def is_ion = het_info ? het_info.is_ion : rec.is_ion

                def pose_meta = [
                    id: p_name,
                    name: p_name,
                    res_name: rec.res_name,
                    inchikey: p_inchikey,
                    charge: rec.charge as int,
                    is_ion: is_ion
                ]
                def sdf_file = json_file.parent.resolve(rec.sdf)
                def mae_file = json_file.parent.resolve(rec.mae)
                [ pose_meta, sdf_file, mae_file ]
            }
        }

    emit:
    prepared_mae     = SPLIT_PREPARED.out.prepared_mae
    membrane_json    = ch_mem_json
    poses_to_param   = ch_poses_to_param
}
