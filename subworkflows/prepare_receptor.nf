/*
 * subworkflows/prepare_receptor.nf
 * Decoupled Receptor preparation pipeline:
 * 1. EXTRACT_COMPLEX (Schrödinger): complex (.cif/.pdb) -> raw_protein.pdb + raw_native.mae + raw_native.sdf
 * 2. If membrane:
 *      IMMERS (Pure PPM 3.0 CLI): raw_protein.pdb -> oriented_protein.pdb + membrane_info
 *      ALIGN_LIGAND (Schrödinger): raw_protein.pdb + oriented_protein.pdb + raw_native.mae -> ref_native.mae + ref_native.sdf
 * 3. PREPWIZARD (Schrödinger): oriented (or raw) protein -> prepared_mae + prepared_pdb
 * 4. If membrane:
 *      APPLY_MEMBRANE (Schrödinger Prime): prepared_mae + membrane_info -> membrane_mae
 */

include { EXTRACT_COMPLEX } from "../modules/schrodinger/extract_complex"
include { IMMERS          } from "../modules/membrane/immers"
include { ALIGN_LIGAND    } from "../modules/schrodinger/align_ligand"
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
    def extract_script = file("${projectDir}/bin/extract_complex.py")
    def parse_script   = file("${projectDir}/bin/parse_immers.py")
    def align_script   = file("${projectDir}/bin/align_ligand.py")
    def apply_script   = file("${projectDir}/bin/apply_membrane.py")

    // 1. Structure extraction: Separate raw protein from native ligand/cofactors
    EXTRACT_COMPLEX(
        ch_input,
        extract_script,
        native_ligand,
        cofactors
    )

    if (membrane_enabled) {
        // 2. Pure PPM 3.0 membrane orientation calculation (standalone, no Schrödinger)
        IMMERS(
            EXTRACT_COMPLEX.out.raw_protein_pdb,
            parse_script,
            membrane_type,
            membrane_topo,
            membrane_chains
        )

        // 3. Align native reference ligand to the oriented receptor frame
        ch_align_in = EXTRACT_COMPLEX.out.raw_protein_pdb
            .join(IMMERS.out.oriented_pdb)
            .join(EXTRACT_COMPLEX.out.raw_native_mae)
        ALIGN_LIGAND(ch_align_in, align_script)

        ch_prep_in          = IMMERS.out.oriented_pdb
        ch_final_ref_native = ALIGN_LIGAND.out.ref_native_mae
        ch_final_ref_sdf    = ALIGN_LIGAND.out.ref_native_sdf
        ch_mem_json         = IMMERS.out.membrane_json
        ch_mem_txt          = IMMERS.out.membrane_txt
    } else {
        ch_prep_in          = EXTRACT_COMPLEX.out.raw_protein_pdb
        ch_final_ref_native = EXTRACT_COMPLEX.out.raw_native_mae
        ch_final_ref_sdf    = EXTRACT_COMPLEX.out.raw_native_sdf
        ch_mem_json         = Channel.empty()
        ch_mem_txt          = Channel.empty()
    }

    // 4. Schrodinger Protein Preparation Wizard
    PREPWIZARD(
        ch_prep_in,
        target_ph,
        ph_tol,
        min_rmsd,
        cap_termini
    )

    ch_final_mae = PREPWIZARD.out.prepared_mae
    ch_final_pdb = PREPWIZARD.out.prepared_pdb

    // 5. Optional implicit membrane property injection for Prime
    if (membrane_enabled) {
        ch_apply_in = PREPWIZARD.out.prepared_mae
            .join(ch_mem_txt)
            .map { meta, mae, txt -> [ meta, mae, txt, apply_script ] }
            
        APPLY_MEMBRANE(ch_apply_in)
        ch_final_mae = APPLY_MEMBRANE.out.membrane_mae
    }

    emit:
    prepared_pdb     = ch_final_pdb
    prepared_mae     = ch_final_mae
    ref_native_mae   = ch_final_ref_native
    ref_native_sdf   = ch_final_ref_sdf
    membrane_json    = ch_mem_json
}
