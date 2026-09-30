/*
 * subworkflows/setup_md_system.nf
 * Full MD System Build: Extraction -> Amberization -> Membrane/Solvation -> GROMACS conversion -> Restraints
 */

include { ANTECHAMBER    } from '../modules/amber/antechamber'
include { MERGE_PDB      } from '../modules/amber/merge_pdb'
include { PACKMOL_MEMGEN } from '../modules/membrane/packmol'
include { TLEAP          } from '../modules/amber/tleap'
include { PARMED         } from '../modules/amber/parmed'
include { DECLASH        } from '../modules/membrane/declash'
include { NDX_POSRE      } from '../modules/gromacs/ndx_posre'

process EXTRACT_MD_SYSTEM {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(rec_file), path(lig_file)
    val(cofactors)

    output:
    tuple val(meta), path("receptor_dry.pdb"), path("organics.json"), path("org_files/*"), emit: extracted

    script:
    def cof_list = cofactors instanceof List ? cofactors.join(' ') : cofactors
    def cof_arg = cof_list ? "--cofactors ${cof_list}" : ""
    def lig_arg = meta.native_ligand ? "--lig_name ${meta.native_ligand}" : ""
    """
    schrodinger run ${projectDir}/bin/extract_md_files.py \\
        --rec_file ${rec_file} \\
        --lig_file ${lig_file} \\
        ${lig_arg} \\
        ${cof_arg} \\
        --out_dir .
    """

    stub:
    """
    touch receptor.pdb
    touch ligand.mol2
    echo "0" > lig_charge.txt
    """
}

workflow SETUP_MD_SYSTEM {
    take:
    ch_targets        // channel: [ val(meta), path(receptor_pdb), path(ligand_file), val(cofactors_list) ]
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
    ch_ext_in = ch_targets.map { meta, rec_pdb, lig_file, cofs -> [ meta, rec_pdb, lig_file ] }
    EXTRACT_MD_SYSTEM(ch_ext_in, ch_targets.map { it[3] })
    
    // 2. Parse organic components for parallel Antechamber
    ch_ante_in = EXTRACT_MD_SYSTEM.out.extracted.flatMap { meta, rec_pdb, json_file, org_files ->
        def org_list = org_files instanceof List ? org_files : [org_files]
        def data = new groovy.json.JsonSlurper().parseText(json_file.text)
        data.collect { org ->
            def mol_meta = meta + [mol_id: org.name]
            def m_file = org_list.find { it.name == org.mol2 }
            def c_file = org_list.find { it.name == org.charge }
            [ mol_meta, m_file, c_file ]
        }
    }

    ANTECHAMBER(ch_ante_in)

    ch_ante_grouped = ANTECHAMBER.out.gaff_params
        .map { mol_meta, mol2, frcmod, pdb ->
            [ [id: mol_meta.id, project: mol_meta.project], mol2, frcmod, pdb ]
        }
        .groupTuple(by: 0, sort: true)

    // 3. Clean Bash Merge of standardized protein and ALL organics
    ch_merge_in = EXTRACT_MD_SYSTEM.out.extracted
        .map { meta, rec_pdb, json_file, dir -> [ meta, rec_pdb ] }
        .join(ch_ante_grouped, remainder: true)
        .map { tuple ->
            def meta = tuple[0]
            def rec_pdb = tuple[1]
            def pdbs = (tuple.size() > 4 && tuple[4] != null) ? tuple[4] : []
            [ meta, rec_pdb, pdbs ]
        }

    MERGE_PDB(ch_merge_in)

    // 4. Membrane Assembly or Solvation
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

    // 5. tLEaP topology building
    ch_tleap_in = ch_tleap_pdb
        .join(ch_ante_grouped, remainder: true)
        .map { tuple ->
            def meta = tuple[0]
            def pdb = tuple[1]
            def mol2s = (tuple.size() > 2 && tuple[2] != null) ? tuple[2] : []
            def frcmods = (tuple.size() > 3 && tuple[3] != null) ? tuple[3] : []
            [ meta, pdb, mol2s, frcmods ]
        }
    
    TLEAP(ch_tleap_in, water_model, is_membrane, wat_pad)

    // 6. ParmEd conversion to GROMACS
    PARMED(TLEAP.out.amber_system)

    // 7. Declash
    if (is_membrane.toString() == 'true') {
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
