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
include { APPLY_MEMBRANE  } from "../modules/membrane/apply_membrane"

workflow PREPARE_RECEPTOR {
    take:
    ch_input          // channel: [ val(meta), path(complex_file), val(frame_idx) ]
    native_ligand     // val: "PKA"
    cofactors         // val: "HEM,MG"
    membrane_enabled  // val: true/false
    membrane_type     // val: "MOM"
    membrane_topo     // val: "in"
    membrane_chains   // val: "A"
    target_ph         // val: 7.4
    ph_tol            // val: 1.0
    min_rmsd          // val: 0.3
    cap_termini       // val: true

    main:
    // 1. Separate models and extract components
    EXTRACT_COMPLEX(
        ch_input,
        native_ligand,
        cofactors,
        membrane_enabled,
        file("${projectDir}/bin/extract_complex.py")
    )

    if (membrane_enabled) {
        // 2. Parallel PPM 3.0 immersion on pure protein
        IMMERS(
            EXTRACT_COMPLEX.out.pure_protein_pdb,
            membrane_type,
            membrane_topo,
            membrane_chains
        )

        // 3. Superimpose and transform entire complex structure (protein + cofactors + ligand)
        ch_orient_in = EXTRACT_COMPLEX.out.pure_protein_pdb
            .join(IMMERS.out.oriented_pdb)
            .join(EXTRACT_COMPLEX.out.complex_mae)

        ORIENT_COMPLEX(
            ch_orient_in,
            native_ligand,
            cofactors,
            file("${projectDir}/bin/orient_complex.py")
        )

        ch_prep_in          = ORIENT_COMPLEX.out.receptor_raw_mae
        ch_final_ref_native = ORIENT_COMPLEX.out.ref_native_mae
        ch_final_ref_sdf    = ORIENT_COMPLEX.out.ref_native_sdf
        ch_mem_json         = IMMERS.out.membrane_json
        ch_mem_txt          = IMMERS.out.membrane_txt
    } else {
        ch_prep_in          = EXTRACT_COMPLEX.out.raw_receptor_mae
        ch_final_ref_native = EXTRACT_COMPLEX.out.ref_native_mae
        ch_final_ref_sdf    = EXTRACT_COMPLEX.out.ref_native_sdf
        ch_mem_json         = Channel.empty()
        ch_mem_txt          = Channel.empty()
    }

    // 4. Protein Preparation Wizard in parallel across Slurm nodes
    PREPWIZARD(
        ch_prep_in,
        target_ph,
        ph_tol,
        min_rmsd,
        cap_termini
    )

    // 5. Implicit membrane property injection for Prime MM-GBSA
    if (membrane_enabled) {
        ch_apply_in = PREPWIZARD.out.prepared_mae
            .join(ch_mem_txt)

        APPLY_MEMBRANE(ch_apply_in, file("${projectDir}/bin/apply_membrane.py"))
        ch_final_mae = APPLY_MEMBRANE.out.prepared_mae
        ch_final_pdb = APPLY_MEMBRANE.out.prepared_pdb
    } else {
        ch_final_mae = PREPWIZARD.out.prepared_mae
        ch_final_pdb = PREPWIZARD.out.prepared_pdb
    }

    emit:
    prepared_pdb     = ch_final_pdb
    prepared_mae     = ch_final_mae
    ref_native_mae   = ch_final_ref_native
    ref_native_sdf   = ch_final_ref_sdf
    membrane_json    = ch_mem_json
}
