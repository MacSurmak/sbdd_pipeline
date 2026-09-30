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

include { COUNT_MODELS    } from "../modules/common/count_models"
include { EXTRACT_COMPLEX } from "../modules/schrodinger/extract_complex"
include { IMMERS          } from "../modules/membrane/immers"
include { ORIENT_COMPLEX  } from "../modules/schrodinger/orient_complex"
include { PREPWIZARD      } from "../modules/schrodinger/prepwizard"
include { APPLY_MEMBRANE   } from "../modules/membrane/apply_membrane"
include { SPLIT_PREPARED   } from "../modules/schrodinger/split_prepared"
include { LIGHT_PREPWIZARD } from "../modules/schrodinger/light_prep"

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
    // 1. Strict branching based on DB prep_mode
    ch_receptor.branch {
        full:  it[0].prep_mode == 'full'
        light: it[0].prep_mode == 'light'
    }.set { routed_targets }

    // ==========================================
    // ROUTE A: FULL PREPARATION (Crystals/Raw structures)
    // ==========================================
    COUNT_MODELS(routed_targets.full)

    ch_extracted_in = routed_targets.full
        .join(COUNT_MODELS.out.count_file)
        .flatMap { meta, rec_file, count_txt ->
            def total = count_txt.text.trim() as int
            def models = (1..total).toList()
            models.collect { idx ->
                def suffix = total > 1 ? "_f${idx}" : ""
                def frame_meta = meta + [
                    id: "${meta.id}${suffix}",
                    frame: idx
                ]
                [ frame_meta, rec_file, idx ]
            }
        }

    EXTRACT_COMPLEX(
        ch_extracted_in,
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
    ch_full_prepared_mae = SPLIT_PREPARED.out.prepared_mae

    // ==========================================
    // ROUTE B: LIGHT PREPARATION (MD frames)
    // ==========================================
    LIGHT_PREPWIZARD(routed_targets.light)
    ch_light_prepared_mae = LIGHT_PREPWIZARD.out.prepared_mae


    // ==========================================
    // EMIT RECEIPT FOR CLI DATABASE
    // ==========================================
    ch_all_maes = ch_full_prepared_mae.mix(ch_light_prepared_mae)

    ch_all_maes.map { meta, mae -> "${meta.name}\t${mae.toAbsolutePath()}" }
               .collectFile(name: 'prep_receipt.tsv', storeDir: params.outdir, newLine: true)

    emit:
    prepared_mae     = ch_all_maes
    ref_native_mae   = SPLIT_PREPARED.out.ref_native_mae
    ref_native_sdf   = SPLIT_PREPARED.out.ref_native_sdf
    membrane_json    = ch_mem_json
}
