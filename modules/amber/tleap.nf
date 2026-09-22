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
    import math

    with open("${complex_pdb}", "r") as f:
        lines = f.readlines()
    cleaned = []
    last_chain = "A"
    cys_atoms = []

    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")) and len(l) >= 22:
            ch = l[21].strip()
            if ch: last_chain = ch
            if l[12:16].strip() == "SG" and l[17:20].strip() in ["CYS", "CYX", "CYM"]:
                cys_atoms.append({
                    "chain": l[21],
                    "resSeq": int(l[22:26]),
                    "x": float(l[30:38]),
                    "y": float(l[38:46]),
                    "z": float(l[46:54])
                })
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

    # Detect disulfide bonds within 2.3 A
    bonds = []
    for i in range(len(cys_atoms)):
        for j in range(i + 1, len(cys_atoms)):
            a, b = cys_atoms[i], cys_atoms[j]
            if a["resSeq"] == b["resSeq"] and a["chain"] == b["chain"]: continue
            d = math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)
            if d <= 2.3:
                bonds.append(f"bond complex.{a['resSeq']}.SG complex.{b['resSeq']}.SG")

    with open("disulfide_bonds.txt", "w") as fb:
        fb.write("\n".join(bonds))
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
EOF

    cat disulfide_bonds.txt >> tleap.in

    cat << EOF >> tleap.in
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

