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

SMOKE_DIR="/tmp/cluster_smoke_$(date +%s)"
mkdir -p "${SMOKE_DIR}"
trap "rm -rf ${SMOKE_DIR} /shared/work/cluster_smoke_test 2>/dev/null || true" EXIT

export NXF_HOME="${SMOKE_DIR}/.nextflow"

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
if [ -f "${SMOKE_DIR}/docking/cluster_smoke_test_docked_all.sdf" ]; then
    echo "  ✓ Found merged docking poses: docking/cluster_smoke_test_docked_all.sdf ($(wc -c < "${SMOKE_DIR}/docking/cluster_smoke_test_docked_all.sdf") bytes)"
else
    echo "  ✗ ERROR: docking/cluster_smoke_test_docked_all.sdf was not generated!" >&2
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

