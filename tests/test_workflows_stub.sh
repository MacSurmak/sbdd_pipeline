#!/usr/bin/env bash
#
# tests/test_workflows_stub.sh
# Automated Nextflow DSL2 DAG & Process Verification on Slurm Cluster in Stub Mode
# Fast verification (~1 minute) by validating multi-frame/multi-replicate DAGs with 2 ensemble models.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
MAIN_NF="${ROOT_DIR}/main.nf"

echo "=========================================================="
echo " Starting Nextflow DSL2 Workflow Stub Tests on Slurm"
echo " Framework Root: ${ROOT_DIR}"
echo "=========================================================="

TMP_OUT="/shared/work/sbdd_stub_tests_$(date +%s)"
mkdir -p "${TMP_OUT}"
trap "rm -rf ${TMP_OUT}" EXIT

# 1. Test PREPARE workflow
echo -e "\n[1/5] Testing PREPARE workflow on Slurm cluster (-stub)..."
nextflow -log "${TMP_OUT}/nextflow_prep.log" run "${MAIN_NF}" \
    -stub \
    -work-dir "${TMP_OUT}/work_prep" \
    --entry prepare \
    --total_models 2 \
    --project_name test_stub_prep \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --library "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/prep" \
    -profile slurm

# 2. Test SCREEN workflow (GNINA)
echo -e "\n[2/5] Testing SCREEN workflow (GNINA) on Slurm cluster (-stub)..."
nextflow -log "${TMP_OUT}/nextflow_screen_gnina.log" run "${MAIN_NF}" \
    -stub \
    -work-dir "${TMP_OUT}/work_screen_gnina" \
    --entry screen \
    --total_models 2 \
    --dock_engine gnina \
    --project_name test_stub_screen_gnina \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --library "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/screen_gnina" \
    -profile slurm

# 3. Test SCREEN workflow (Both: Glide + GNINA Rescore + MM-GBSA 2D)
echo -e "\n[3/5] Testing SCREEN workflow (Both: Glide + GNINA + 2D MM-GBSA) on Slurm (-stub)..."
nextflow -log "${TMP_OUT}/nextflow_screen_both.log" run "${MAIN_NF}" \
    -stub \
    -work-dir "${TMP_OUT}/work_screen_both" \
    --entry screen \
    --total_models 2 \
    --dock_engine both \
    --project_name test_stub_screen_both \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --library "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/screen_both" \
    -profile slurm

# 4. Test DYNAMICS workflow
echo -e "\n[4/5] Testing DYNAMICS workflow on Slurm cluster (-stub)..."
nextflow -log "${TMP_OUT}/nextflow_dyn.log" run "${MAIN_NF}" \
    -stub \
    -work-dir "${TMP_OUT}/work_dyn" \
    --entry dynamics \
    --md_replicates 1 \
    --project_name test_stub_dyn \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --docked_poses "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/dynamics" \
    -profile slurm

# 5. Test SBDD master end-to-end workflow
echo -e "\n[5/5] Testing Master SBDD workflow on Slurm cluster (-stub)..."
nextflow -log "${TMP_OUT}/nextflow_all.log" run "${MAIN_NF}" \
    -stub \
    -work-dir "${TMP_OUT}/work_all" \
    --entry all \
    --total_models 2 \
    --md_replicates 1 \
    --project_name test_stub_all \
    --receptor "${SCRIPT_DIR}/../tests/fixtures/dummy_rec.pdb" \
    --library "${SCRIPT_DIR}/../tests/fixtures/dummy_lib.sdf" \
    --outdir "${TMP_OUT}/all" \
    -profile slurm

echo -e "\n=========================================================="
echo " ✓ ALL 5 NEXTFLOW WORKFLOW DAGs PASSED ON SLURM IN STUB MODE!"
echo "=========================================================="
