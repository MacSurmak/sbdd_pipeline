/*
 * modules/amber/tleap.nf
 * Amber tLEaP System Topology and Coordinate Builder
 */

process TLEAP {
    tag "${meta.id}"
    label 'process_single'

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
    # 1. Standardize Schrödinger terminal caps and protonation to AMBER format using ParmEd API
    amber_cap_fixer.py ${complex_pdb} complex_clean.pdb

    # 2. Map SS-links dynamically
    python3 -c "
import sys
from pathlib import Path
# Append framework bin to path so we can import topology_prep
sys.path.append('${projectDir}/bin/utils')
try:
    # Try importing from the bin directory structure if copied
    sys.path.append(str(Path('${projectDir}/bin')))
except Exception:
    pass

def generate_bonds(pdb):
    import math
    cys = []
    with open(pdb) as f:
        for l in f:
            if l.startswith(('ATOM','HETATM')) and l[12:16].strip() == 'SG' and l[17:20].strip() in ['CYS','CYX']:
                cys.append({'c': l[21], 'r': int(l[22:26]), 'x': float(l[30:38]), 'y': float(l[38:46]), 'z': float(l[46:54])})
    bonds = []
    for i in range(len(cys)):
        for j in range(i+1, len(cys)):
            a, b = cys[i], cys[j]
            if a['r'] == b['r'] and a['c'] == b['c']: continue
            if math.sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2) <= 2.3:
                bonds.append(f'bond complex.{a[\\"r\\"]}.SG complex.{b[\\"r\\"]}.SG')
    return '\\n'.join(bonds)

open('apply_bonds.in', 'w').write(generate_bonds('complex_clean.pdb') + '\\n')
"

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
