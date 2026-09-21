/*
 * modules/amber/tleap.nf
 * Amber tLEaP System Topology and Coordinate Builder
 */

process TLEAP {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(complex_pdb), path(unk_mol2), path(unk_frcmod)
    val(water_model)
    val(is_membrane)
    val(box_padding)

    output:
    tuple val(meta), path("complex.prmtop"), path("complex.inpcrd"), emit: amber_system
    tuple val(meta), path("complex_tleap.pdb"), emit: tleap_pdb

    script:
    def solvate_cmd = !is_membrane ? "solvateBox complex OPCBOX ${box_padding}" : ""
    """
    # Clean artificial TER around terminal caps
    python3 - << 'EOF'
    with open("${complex_pdb}", "r") as f:
        lines = f.readlines()
    cleaned = []
    last_chain = "A"
    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")) and len(l) >= 22:
            ch = l[21].strip()
            if ch: last_chain = ch
        if l.startswith("TER"):
            if i > 0 and "ACE" in lines[i-1]: continue
            if i + 1 < len(lines) and "NME" in lines[i+1]: continue
        if ("ACE" in l or "NME" in l) and l.startswith(("ATOM", "HETATM")):
            l = l[:21] + last_chain + l[22:]
        if "ASP" in l and ("HD2" in l or "2HD" in l): l = l[:17] + "ASH" + l[20:]
        elif "GLU" in l and ("HE2" in l or "2HE" in l): l = l[:17] + "GLH" + l[20:]
        cleaned.append(l)
    with open("complex_clean.pdb", "w") as f:
        f.writelines(cleaned)
    EOF

    cat << EOF > tleap.in
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.lipid21
source leaprc.water.opc
loadAmberParams frcmod.ionslm_126_opc
UNK = loadmol2 ${unk_mol2}
loadAmberParams ${unk_frcmod}
complex = loadPdb complex_clean.pdb
set default PBRadii mbondi3
addPdbResMap { {"Zn2+" "ZN"} {"Ca2+" "CA"} {"Mg2+" "MG"} {"Na+" "NA"} {"K+" "K"} {"Cl-" "CL"} }
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

