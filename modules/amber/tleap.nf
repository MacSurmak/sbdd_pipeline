/*
 * modules/amber/tleap.nf
 * Amber tLEaP System Topology and Coordinate Builder
 */

process TLEAP {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(complex_pdb), path(mol2s), path(frcmods)
    val(water_model)
    val(is_membrane)
    val(box_padding)

    output:
    tuple val(meta), path("complex.prmtop"), path("complex.inpcrd"), emit: amber_system
    tuple val(meta), path("complex_tleap.pdb"), emit: tleap_pdb

    script:
    def is_mem = (is_membrane.toString() == 'true')
    def is_tip3p = (water_model.toString().toLowerCase() == 'tip3p')
    def water_leaprc = is_tip3p ? "source leaprc.water.tip3p" : "source leaprc.water.opc"
    def ions_params = is_tip3p ? "loadAmberParams frcmod.ionslm_126_tip3p" : "loadAmberParams frcmod.ionslm_126_opc"
    def box_name = is_tip3p ? "TIP3PBOX" : "OPCBOX"
    def solvate_cmd = !is_mem ? "solvateBox complex ${box_name} ${box_padding}" : ""
    """
    python3 ${projectDir}/bin/amber_cap_fixer.py ${complex_pdb} complex_clean.pdb --bonds_out apply_bonds.in

    cat << EOF > tleap.in
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.lipid21
${water_leaprc}
${ions_params}
EOF

    # Dynamically load all organics
    for m in ${mol2s}; do
        name=\$(basename "\$m" _gaff2.mol2)
        echo "\$name = loadmol2 \$m" >> tleap.in
    done
    for f in ${frcmods}; do
        echo "loadAmberParams \$f" >> tleap.in
    done

    cat << EOF >> tleap.in
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
