#!/usr/bin/env bash
#
# tests/cluster_smoke/run_smoke_test.sh
# End-to-End Direct MD Smoke Test on Slurm HPC
# Validates initialization, structural prep, protocol generation, and MD on a raw PDB.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "=========================================================="
echo " Starting End-to-End Direct MD Smoke Test"
echo " Framework Root: ${ROOT_DIR}"
echo "=========================================================="

SMOKE_DIR="/shared/smoke_md_$(date +%s)"
mkdir -p "${SMOKE_DIR}"

cleanup() {
    EXIT_STATUS=$?
    if [ $EXIT_STATUS -ne 0 ]; then
        echo -e "\n[DEBUG] Failure detected. Preserving workdir for inspection: ${SMOKE_DIR}"
    else
        echo "  ✓ Cleaning up ${SMOKE_DIR}"
        rm -rf "${SMOKE_DIR}"
    fi
}
trap cleanup EXIT

export NXF_HOME="${SMOKE_DIR}/.nextflow"

echo "[1/4] Preparing data and configurations..."
cd "${SMOKE_DIR}"

# Initialize project structure
sbdd init . > /dev/null

# Disable membrane for fast testing
sed -i 's/enabled: true/enabled: false/g' configs/membrane.yaml

# Download 1STP and extract exactly one frame
python3 -c '
import urllib.request
urllib.request.urlretrieve("https://files.rcsb.org/download/1STP.pdb", "1stp_raw.pdb")
with open("1stp_raw.pdb", "r") as f:
    lines = [l for l in f if l.startswith(("ATOM", "HETATM", "TER"))]
with open("1stp.pdb", "w") as f:
    f.writelines(lines)
'

# Create realistic manifest
cat <<EOF > configs/manifest.yaml
targets:
  - id: "1stp_test"
    structure: "1stp.pdb"
    native_ligand: "BTN"
EOF

# Create simplified MD config (no membrane, tiny times)
cat <<EOF > configs/md.yaml
md_simulation:
  enabled: true
  pbc_center: "Protein"
  pbc_extract: "COMPLEX"
  restrain_groups: ["Protein", "UNK"]
  system_build:
    water_model: "opc"
    box_padding: 8.0
    neutralize: true
    salt_conc: 0.15
  production:
    total_time_ns: 0.002
    chunk_size_ns: 0.002
    replicates: 1
EOF

echo -e "\n[2/4] Running 'sbdd prep' (Standardization & Extraction)..."
sbdd prep

if [ ! -f "prepared/1stp_test_prepared.pdb" ]; then
    echo "✗ ERROR: prepared/1stp_test_prepared.pdb not found!"
    exit 1
fi
if [ ! -f "prepared/1stp_test_ref_native.sdf" ]; then
    echo "✗ ERROR: prepared/1stp_test_ref_native.sdf not found!"
    exit 1
fi
echo "  ✓ Structure properly prepared and native ligand extracted."

echo -e "\n[3/4] Initializing MD protocol and optimizing for speed..."
sbdd md-init > /dev/null

# Modify protocol for a lightning fast smoke test (50 steps per stage)
sed -i 's/nsteps = 50000/nsteps = 50/g' configs/md_protocol/*.mdp
sed -i 's/nsteps = 100000/nsteps = 50/g' configs/md_protocol/*.mdp
sed -i 's/nsteps = 500000/nsteps = 50/g' configs/md_protocol/*.mdp

# Edit the playlist to only run Minimization -> NVT -> Production
echo "minim.mdp" > configs/md_protocol/protocol.txt
echo "nvt.mdp" >> configs/md_protocol/protocol.txt
echo "  ✓ Playlist (protocol.txt) successfully generated and customized."

echo -e "\n[4/4] Running 'sbdd md' using generated protocol..."
sbdd md

echo -e "\n[Verifying MD Outputs]..."
MD_XTC=$(find "md" -name "md_1.xtc" | head -n 1)
if [ -n "$MD_XTC" ] && [ -f "$MD_XTC" ]; then
    echo "  ✓ Found Production MD trajectory: $MD_XTC"
else
    echo "  ✗ ERROR: MD trajectory was not generated!"
    exit 1
fi

SYS_TPR=$(find "md" -name "sys_ref.tpr" | head -n 1)
if [ -n "$SYS_TPR" ] && [ -f "$SYS_TPR" ]; then
    echo "  ✓ Found system topology: $SYS_TPR"
else
    echo "  ✗ ERROR: sys_ref.tpr was not generated!"
    exit 1
fi

echo -e "\n=========================================================="
echo " ✓ DIRECT MD SMOKE TEST PASSED!"
echo "=========================================================="
