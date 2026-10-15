/*
 * subworkflows/setup_md_system.nf
 * Full MD System Build: Extraction -> Amberization -> Membrane/Solvation -> GROMACS conversion -> Restraints
 */

include { MERGE_PDB         } from '../modules/amber/merge_pdb'
include { PACKMOL_MEMGEN    } from '../modules/membrane/packmol'
include { TLEAP             } from '../modules/amber/tleap'
include { PARMED            } from '../modules/amber/parmed'
include { DECLASH           } from '../modules/membrane/declash'
include { NDX_POSRE         } from '../modules/gromacs/ndx_posre'

workflow SETUP_MD_SYSTEM {
    take:
    ch_targets        // channel: [ meta, rec_pdb, lig_sdf, extra_sdfs, mol2s, frcmods ]
    is_membrane       // val: true/false
    upper_lipids      // val: "POPC:POPE:CHL1"
    upper_ratio       // val: "5:4:1"
    lower_lipids      // val: "POPC:POPE:POPS:CHL1"
    lower_ratio       // val: "4:4:1:1"
    lipid_pad         // val: 10.0
    wat_pad           // val: 15.0
    salt_conc         // val: 0.15
    water_leaprc
    water_box
    box_type
    ss_cutoff
    water_names
    ion_names
    lipid_names
    restraints_fc
    declash_cutoff
    declash_niter
    declash_damping
    declash_max_move
    declash_jitter
    declash_tol
    protein_ff
    gaff_ff
    lipid_ff

    main:
    // 1. Clean merge of standardized apo protein with active ligand and cofactors/ions
    ch_merge_in = ch_targets.map { meta, rec_pdb, lig_sdf, extra_sdfs, mol2s, frcmods ->
        [ meta, rec_pdb, lig_sdf, extra_sdfs ]
    }

    MERGE_PDB(ch_merge_in)

    // 2. Membrane Assembly or Solvation
    if (is_membrane.toString() == 'true') {
        PACKMOL_MEMGEN(
            MERGE_PDB.out.complex_pdb,
            upper_lipids,
            upper_ratio,
            lower_lipids,
            lower_ratio,
            lipid_pad,
            wat_pad,
            salt_conc
        )
        ch_tleap_pdb = PACKMOL_MEMGEN.out.mem_complex_pdb
    } else {
        ch_tleap_pdb = MERGE_PDB.out.complex_pdb
    }

    // 3. tLEaP topology building
    ch_tleap_in = ch_tleap_pdb
        .join(ch_targets.map { meta, rec_pdb, lig_sdf, extra_sdfs, mol2s, frcmods -> [meta, mol2s, frcmods] })
        .map { meta, pdb, mol2s, frcmods ->
            [ meta, pdb, mol2s, frcmods ]
        }
    
    TLEAP(
        ch_tleap_in,
        water_leaprc,
        water_box,
        box_type,
        is_membrane,
        wat_pad,
        ss_cutoff,
        water_names,
        salt_conc,
        protein_ff,
        gaff_ff,
        lipid_ff
    )

    // 6. ParmEd conversion to GROMACS
    PARMED(TLEAP.out.amber_system)

    // 7. Declash (Molecules-based)
    if (is_membrane.toString() == 'true') {
        ch_declash_in = PARMED.out.gmx_system.map { meta, topol, gro -> [ meta, gro, topol ] }
        DECLASH(
            ch_declash_in,
            declash_cutoff,
            declash_niter,
            declash_damping,
            declash_max_move,
            declash_jitter,
            declash_tol
        )
        ch_final_gro = DECLASH.out.declashed_gro
    } else {
        ch_final_gro = PARMED.out.gmx_system.map { meta, topol, gro -> [ meta, gro ] }
    }

    // 8. Index & Position Restraints
    ch_ndx_in = ch_final_gro.join(PARMED.out.gmx_system.map { meta, topol, gro -> [ meta, topol ] })
    NDX_POSRE(ch_ndx_in, restraints_fc, water_names, ion_names, lipid_names)

    ch_system_ready = ch_final_gro.join(NDX_POSRE.out.restraints)
        .map { meta, gro, ndx, topol, posre_files ->
            [ meta, gro, topol, ndx, posre_files ]
        }

    emit:
    system_ready = ch_system_ready
}
