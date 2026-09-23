/*
 * subworkflows/setup_md_system.nf
 * Full MD System Build: Extraction -> Amberization -> Membrane/Solvation -> GROMACS conversion -> Restraints
 */

include { ANTECHAMBER    } from '../modules/amber/antechamber'
include { PACKMOL_MEMGEN } from '../modules/membrane/packmol'
include { TLEAP          } from '../modules/amber/tleap'
include { PARMED         } from '../modules/amber/parmed'
include { DECLASH        } from '../modules/membrane/declash'
include { NDX_POSRE      } from '../modules/gromacs/ndx_posre'

process EXTRACT_MD_SYSTEM {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(lig_file), val(lig_title), path(rec_file)

    output:
    tuple val(meta), path("receptor.pdb"), path("ligand.mol2"), path("lig_charge.txt"), emit: extracted

    script:
    """
    extract_md_files.py \\
        --lig_file ${lig_file} \\
        --lig_title "${lig_title}" \\
        --rec_file ${rec_file} \\
        --out_dir .
    """

    stub:
    """
    touch receptor.pdb
    touch ligand.mol2
    echo "0" > lig_charge.txt
    """
}

process MERGE_AMBER_COMPLEX {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_pdb), path(unk_mol2)

    output:
    tuple val(meta), path("complex.pdb"), emit: complex_pdb

    script:
    """
    amber_fix_pdb.py ${rec_pdb} fixed_receptor.pdb
    merge_complex.py fixed_receptor.pdb ${unk_mol2} complex.pdb
    """

    stub:
    """
    touch complex.pdb
    """
}

workflow SETUP_MD_SYSTEM {
    take:
    ch_targets        // channel: [ val(meta), path(docked_pose), val(hit_title), path(receptor_file) ]
    is_membrane       // val: true/false
    upper_lipids      // val: "POPC:POPE:CHL1"
    upper_ratio       // val: "5:4:1"
    lower_lipids      // val: "POPC:POPE:POPS:CHL1"
    lower_ratio       // val: "4:4:1:1"
    lipid_pad         // val: 10.0
    wat_pad           // val: 15.0
    salt_conc         // val: 0.15
    water_model       // val: "opc"

    main:
    // 1. Extract files
    EXTRACT_MD_SYSTEM(ch_targets)
    
    ch_ante_in = EXTRACT_MD_SYSTEM.out.extracted.map { meta, rec_pdb, lig_mol2, charge_txt ->
        [ meta, lig_mol2, charge_txt ]
    }

    // 2. Antechamber GAFF2 / AM1-BCC
    ANTECHAMBER(ch_ante_in)

    // 3. Merge Amber Complex: join receptor.pdb with UNK_gaff2.mol2
    ch_merge_in = EXTRACT_MD_SYSTEM.out.extracted
        .map { meta, rec_pdb, lig_mol2, charge_txt -> [ meta, rec_pdb ] }
        .join(ANTECHAMBER.out.gaff_params)
        .map { meta, rec_pdb, unk_mol2, unk_frcmod -> [ meta, rec_pdb, unk_mol2 ] }

    MERGE_AMBER_COMPLEX(ch_merge_in)

    // 4. Membrane Assembly or Solvation
    if (is_membrane) {
        PACKMOL_MEMGEN(
            MERGE_AMBER_COMPLEX.out.complex_pdb,
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
        ch_tleap_pdb = MERGE_AMBER_COMPLEX.out.complex_pdb
    }

    // 5. tLEaP topology building
    ch_tleap_in = ch_tleap_pdb.join(ANTECHAMBER.out.gaff_params)
    TLEAP(ch_tleap_in, water_model, is_membrane, wat_pad)

    // 6. ParmEd conversion to GROMACS
    PARMED(TLEAP.out.amber_system)

    // 7. Declashing if membrane system
    if (is_membrane) {
        DECLASH(PARMED.out.gmx_system)
        ch_final_gro = DECLASH.out.declashed_gro
    } else {
        ch_final_gro = PARMED.out.gmx_system.map { meta, topol, gro -> [ meta, gro ] }
    }

    // 8. Index & Position Restraints
    ch_ndx_in = ch_final_gro.join(PARMED.out.gmx_system.map { meta, topol, gro -> [ meta, topol ] })
    NDX_POSRE(ch_ndx_in)

    ch_system_ready = ch_final_gro.join(NDX_POSRE.out.restraints)
        .map { meta, gro, ndx, topol, posre_files ->
            [ meta, gro, topol, ndx, posre_files ]
        }

    emit:
    system_ready = ch_system_ready
}

