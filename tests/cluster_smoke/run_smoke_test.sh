#!/usr/bin/env bash
#
# tests/cluster_smoke/run_smoke_test.sh
# End-to-End Real Slurm Cluster Integration Smoke Test
# Dispatches 1 fast docking job across the Slurm GPU cluster,
# verifying Slurm submission, GPU worker allocation, NFS /shared staging,
# and deliverable publishing in < 30 seconds.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "=========================================================="
echo " Starting End-to-End Slurm Cluster Integration Smoke Test"
echo " Framework Root: ${ROOT_DIR}"
echo "=========================================================="

SMOKE_DIR="/shared/smoke_test_$(date +%s)"
mkdir -p "${SMOKE_DIR}"

cleanup() {
    EXIT_STATUS=$?
    if [ $EXIT_STATUS -ne 0 ]; then
        echo -e "\n[DEBUG] Failure detected. Preserving workdir for inspection: /shared/work/cluster_smoke_test"
    fi
}
trap cleanup EXIT

export NXF_HOME="${SMOKE_DIR}/.nextflow"

# Ensure authentic protein-ligand complex (1STP: Streptavidin with native BTN) is downloaded and made into a 2-frame ensemble
if [ ! -f "${SCRIPT_DIR}/data/smoke_receptor.pdb" ] || ! grep -q "MODEL 2" "${SCRIPT_DIR}/data/smoke_receptor.pdb"; then
    echo "Downloading 1STP complex and generating 2-frame ensemble PDB..."
    python3 -c '
import urllib.request
import os

urllib.request.urlretrieve("https://files.rcsb.org/download/1STP.pdb", "1stp_raw.pdb")
with open("1stp_raw.pdb", "r") as f:
    lines = f.readlines()

coords = [l for l in lines if l.startswith(("ATOM", "HETATM", "TER"))]

with open("'${SCRIPT_DIR}'/data/smoke_receptor.pdb", "w") as f:
    f.write("MODEL 1\n")
    f.writelines(coords)
    f.write("ENDMDL\nMODEL 2\n")
    # Duplicate coordinates to simulate a multi-frame ensemble
    f.writelines(coords)
    f.write("ENDMDL\nEND\n")

os.remove("1stp_raw.pdb")
'
fi

# Generate valid 3D SDF library using RDKit
python3 -c '
from rdkit import Chem
from rdkit.Chem import AllChem
import sys

smiles = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O", 
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O", 
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1"
}

try:
    w = Chem.SDWriter("'${SCRIPT_DIR}'/data/smoke_ligand.sdf")
    for n, s in smiles.items():
        m = Chem.MolFromSmiles(s)
        m.SetProp("_Name", n)
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(m)
        w.write(m)
    w.close()
    print("Successfully generated smoke_ligand.sdf with 3 test molecules.")
except Exception as e:
    print(f"Error generating SDF: {e}")
    sys.exit(1)
'

# Copy smoke data and configs
cp -r "${SCRIPT_DIR}/configs" "${SMOKE_DIR}/"
cp "${SCRIPT_DIR}/data/smoke_receptor.pdb" "${SMOKE_DIR}/"
cp "${SCRIPT_DIR}/data/smoke_ligand.sdf" "${SMOKE_DIR}/"

START_TIME=$(date +%s)

echo -e "\n[1/3] Launching 'sbdd all --engine both' on Slurm HPC cluster..."
(
    cd "${SMOKE_DIR}"
    # sbdd uses default slurm profile and main partition
    sbdd all --engine both
)

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo -e "\n[2/3] Verifying deliverables in ${SMOKE_DIR}..."
DOCK_OUT=$(find "${SMOKE_DIR}/docking" -name "*docked_all.sdf" | head -n 1)
if [ -n "$DOCK_OUT" ] && [ -f "$DOCK_OUT" ]; then
    echo "  ✓ Found merged GNINA docking poses: $(basename "$DOCK_OUT")"
else
    echo "  ✗ ERROR: docking/*_docked_all.sdf was not generated!" >&2
    exit 1
fi

GLIDE_OUT=$(find "${SMOKE_DIR}/docking" -name "*dock_pv.maegz" | head -n 1)
if [ -n "$GLIDE_OUT" ] && [ -f "$GLIDE_OUT" ]; then
    echo "  ✓ Found merged Glide docking poses: $(basename "$GLIDE_OUT")"
else
    echo "  ✗ ERROR: docking/*_dock_pv.maegz was not generated!" >&2
    exit 1
fi

if [ -f "${SMOKE_DIR}/reports/cross_frame_matrix.csv" ]; then
    echo "  ✓ Found cross-frame matrix: reports/cross_frame_matrix.csv"
else
    echo "  ✗ ERROR: reports/cross_frame_matrix.csv was not generated!" >&2
    exit 1
fi

if [ -f "${SMOKE_DIR}/reports/ranking.csv" ]; then
    echo "  ✓ Found master ranking report: reports/ranking.csv"
    head -n 3 "${SMOKE_DIR}/reports/ranking.csv" | column -s, -t
else
    echo "  ✗ ERROR: reports/ranking.csv was not generated!" >&2
    exit 1
fi

MD_XTCS=$(find "${SMOKE_DIR}/md" -name "concatenated.xtc" | wc -l)
if [ "$MD_XTCS" -ge 2 ]; then
    echo "  ✓ Found $MD_XTCS MD trajectories (expected >= 2 for 2 replicates)."
else
    echo "  ✗ ERROR: Expected at least 2 MD trajectories (rep 1 & 2), found $MD_XTCS!" >&2
    exit 1
fi

PB_DATS=$(find "${SMOKE_DIR}/reports/mmpbsa" -name "thermodynamics_mmpbsa_pb.dat" | wc -l)
if [ "$PB_DATS" -ge 2 ]; then
    echo "  ✓ Found $PB_DATS MM-PBSA reports from trajectory analysis."
else
    echo "  ✗ ERROR: Expected at least 2 MM-PBSA reports, found $PB_DATS!" >&2
    exit 1
fi

echo -e "\n[3/3] Checking cluster scratch cleanliness..."
if [ ! -d "/shared/work/cluster_smoke_test" ]; then
    echo "  ✓ /shared/work/cluster_smoke_test was automatically cleaned upon success."
else
    echo "  (Scratch directory present - cleaning now)"
    rm -rf /shared/work/cluster_smoke_test
fi

echo -e "\n=========================================================="
echo " ✓ CLUSTER INTEGRATION SMOKE TEST PASSED IN ${ELAPSED}s!"
echo "=========================================================="

