/*
 * modules/amber/tleap.nf
 * Amber tLEaP System Topology and Coordinate Builder
 */

process TLEAP {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(complex_pdb), path(sslink), path(unk_mol2), path(unk_frcmod)
    val(water_model)
    val(is_membrane)
    val(box_padding)

    output:
    tuple val(meta), path("complex.prmtop"), path("complex.inpcrd"), emit: amber_system
    tuple val(meta), path("complex_tleap.pdb"), emit: tleap_pdb

    script:
    def solvate_cmd = !is_membrane ? "solvateBox complex OPCBOX ${box_padding}" : ""
    """
    # 1. Use official AmberTools ParmEd CLI to rigorously translate Schrodinger cap nomenclature 
    # to Amber standard nomenclature without any manual/custom python string manipulation.
    cat << EOF > fix_caps.in
change resname :NMA NME
change name *@HN H
change name *@1H H
change name :ACE@H H1
change name :ACE@2H H2
change name :ACE@3H H3
change name :NME@H H1
change name :NME@2H H2
change name :NME@3H H3
change name :NME@CA C
change name :NME@CH3 C
change name :NME@1HA H1
change name :NME@2HA H2
change name :NME@3HA H3
outPdb complex_clean.pdb
EOF
    parmed -O -i fix_caps.in ${complex_pdb}

    # 2. Map SS-links from pdb4amber
    if [ -s "${sslink}" ]; then
        sed 's/rec_amber/complex/g' ${sslink} > apply_bonds.in
    else
        touch apply_bonds.in
    fi

    # 3. Build tLEaP topology
    cat << EOF > tleap.in
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.lipid21
source leaprc.water.opc
loadAmberParams frcmod.ionslm_126_opc

# Load ligand definitions first so the PDB loads without warnings
UNK = loadmol2 ${unk_mol2}
loadAmberParams ${unk_frcmod}

complex = loadPdb complex_clean.pdb
set default PBRadii mbondi3
addPdbResMap { {"Zn2+" "ZN"} {"Ca2+" "CA"} {"Mg2+" "MG"} {"Na+" "NA"} {"K+" "K"} {"Cl-" "CL"} }

source apply_bonds.in

${solvate_cmd}
addIonsRand complex Cl- 0
addIonsRand complex K+ 0
charge complex
saveAmberParm complex complex.prmtop complex.inpcrd
savepdb complex complex_tleap.pdb
quit
EOF

    tleap -f tleap.in
    """

    stub:
    """
    touch complex.prmtop
    touch complex.inpcrd
    touch complex_tleap.pdb
    """
}
