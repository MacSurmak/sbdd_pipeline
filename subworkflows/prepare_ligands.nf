/*
 * subworkflows/prepare_ligands.nf
 * Splits SDF library into parallel chunks and/or runs Schrodinger LigPrep
 */

include { LIGPREP } from '../modules/schrodinger/ligprep'

workflow PREPARE_LIGANDS {
    take:
    ch_lig_batches // channel: [ val(meta), path(spec_json), list(path(input_sdfs)) ]
    target_ph      // val: 7.4
    ph_tol         // val: 2.0
    ff_version     // val: 16
    max_states     // val: 32
    outdir         // val: canonical project outdir

    main:
    LIGPREP(ch_lig_batches, target_ph, ph_tol, ff_version, max_states, outdir)

    emit:
    receipt = LIGPREP.out.receipt
}
