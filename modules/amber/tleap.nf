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
    # 1. Standardize Schrödinger terminal caps to AMBER format using ParmEd API
    amber_cap_fixer.py complex_dry.pdb complex_clean.pdb

    # 2. Map SS-links from pdb4amber (update prefix to match our tLEaP variable)
    if [ -s "${sslink}" ]; then
        sed 's/rec_amber/complex/g' ${sslink} > apply_bonds.in
    else
        touch apply_bonds.in
    fi

    # 3. Build tLEaP topology
    cat << EOF > tleap.insource leaprc.protein.ff19SB
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
