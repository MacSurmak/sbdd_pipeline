#!/usr/bin/env bash
#
# tests/test_workflows_stub.sh
# Automated Nextflow DSL2 DAG & Process Verification in Stub Mode
# Runs each standalone workflow and the master pipeline in seconds.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MAIN_NF="${ROOT_DIR}/main.nf"

echo "=========================================================="
echo " Starting Automated Nextflow DSL2 Workflow Stub Tests"
echo " Framework Root: ${ROOT_DIR}"
echo "=========================================================="

TMP_OUT="/tmp/sbdd_stub_tests_$(date +%s)"
mkdir -p "${TMP_OUT}"
trap "rm -rf ${TMP_OUT}" EXIT

# 1. Test PREPARE workflow
echo -e "\n[1/4] Testing PREPARE workflow DAG (-stub)..."
nextflow run "${MAIN_NF}" \
    -stub \
    --entry prepare \
    --project_name test_stub_prep \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --library "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/prep" \
    -profile local

# 2. Test SCREEN workflow
echo -e "\n[2/4] Testing SCREEN workflow DAG (-stub)..."
nextflow run "${MAIN_NF}" \
    -stub \
    --entry screen \
    --project_name test_stub_screen \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --library "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/screen" \
    -profile local

# 3. Test DYNAMICS workflow
echo -e "\n[3/4] Testing DYNAMICS workflow DAG (-stub)..."
nextflow run "${MAIN_NF}" \
    -stub \
    --entry dynamics \
    --project_name test_stub_dyn \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --docked_poses "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/dynamics" \
    -profile local

# 4. Test SBDD master end-to-end workflow
echo -e "\n[4/4] Testing Master SBDD workflow DAG (-stub)..."
nextflow run "${MAIN_NF}" \
    -stub \
    --entry all \
    --project_name test_stub_all \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --library "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/all" \
    -profile local

echo -e "\n=========================================================="
echo " ✓ ALL 4 NEXTFLOW WORKFLOW DAGs PASSED IN STUB MODE!"
echo "=========================================================="

