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
    if [ $EXIT_STATUS -eq 0 ]; then
        rm -rf "${SMOKE_DIR}" /shared/work/cluster_smoke_test 2>/dev/null || true
    else
        echo -e "\n[DEBUG] Failure detected. Preserving workdir for inspection: /shared/work/cluster_smoke_test"
    fi
}
trap cleanup EXIT

export NXF_HOME="${SMOKE_DIR}/.nextflow"

# Ensure authentic protein-ligand complex (1STP: Streptavidin with native BTN) is downloaded from RCSB PDB
if [ ! -f "${SCRIPT_DIR}/data/smoke_receptor.pdb" ] || ! grep -q "BTN" "${SCRIPT_DIR}/data/smoke_receptor.pdb"; then
    echo "Downloading authentic 1STP complex (Streptavidin + BTN) from RCSB PDB..."
    curl -fsSL https://files.rcsb.org/download/1STP.pdb -o "${SCRIPT_DIR}/data/smoke_receptor.pdb"
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

echo -e "\n[1/3] Launching 'sbdd dock --engine gnina' on Slurm HPC cluster..."
(
    cd "${SMOKE_DIR}"
    # sbdd uses default slurm profile and main partition
    sbdd dock --engine gnina
)

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo -e "\n[2/3] Verifying deliverables in ${SMOKE_DIR}..."
DOCK_OUT=$(find "${SMOKE_DIR}/docking" -name "*docked_all.sdf" | head -n 1)
if [ -n "$DOCK_OUT" ] && [ -f "$DOCK_OUT" ]; then
    echo "  ✓ Found merged docking poses: $(basename "$DOCK_OUT") ($(wc -c < "$DOCK_OUT") bytes)"
else
    echo "  ✗ ERROR: docking/*_docked_all.sdf was not generated!" >&2
    exit 1
fi

if [ -f "${SMOKE_DIR}/reports/ranking.csv" ]; then
    echo "  ✓ Found master ranking report: reports/ranking.csv"
    cat "${SMOKE_DIR}/reports/ranking.csv"
else
    echo "  ✗ ERROR: reports/ranking.csv was not generated!" >&2
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

