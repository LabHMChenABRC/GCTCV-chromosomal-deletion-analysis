#!/bin/bash

# --- SLURM Configuration for ref name: DHv4 ---
# NOTE: The --array value is hardcoded based on the known number of samples and intervals.
# If DH.sample.list or the scatter count in 2.mk.interval.sbatch changes,
# this value MUST be updated. A sanity check is performed in the script body.
# written by BoHan Hou
#
# Expected Samples for DHv4: 6
# Expected Intervals: 30
# Total Jobs = 6 * 30 = 180. Array is 0-179.
#
# Before submitting, ensure the log directory exists:
# mkdir -p log/HC.gVCF/DH/DHv4
#
#SBATCH --job-name=gHC.DHv4
#SBATCH --output=log/HC.gVCF/DH/DHv4/%x_%A_%a.log
#SBATCH --error=log/HC.gVCF/DH/DHv4/%x_%A_%a.log
#SBATCH --array=0-179%30
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --partition=intel-g4-al9_short_serial

set -euo pipefail

# --- Configuration ---
PROJECT_NAME="DH"
REF_NAME="DHv4"
SAMPLE_LIST="DH.sample.list"
PLOIDY=3 # Default ploidy
GENERIC_SCRIPT_PATH="./gHC.one.sh" # Assuming it's in the same directory

PROJECT_DIR="/ceph/sharedfs/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/DH"
REF_GENOME="/ceph/work/abrchmc/Reference/mac/Musa_acuminata/DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta.gz"
OUTPUT_BASE_DIR="${PROJECT_DIR}/calling/gVCF/${REF_NAME}"

CORE=${SLURM_CPUS_PER_TASK:-1}

# --- Sanity Check ---
# This section validates that the hardcoded --array value matches the actual file counts.
echo "--- Performing Sanity Check ---"

# Read the list of CRAM files for the specified reference into an array.
# Using mapfile is safer than SAMPLES=($(...)) as it handles paths with spaces.
mapfile -t SAMPLES < <(grep -E "/${REF_NAME}/[^/]+$" "$SAMPLE_LIST")
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "Error: Failed to find any samples for ref name '${REF_NAME}' in ${SAMPLE_LIST}" >&2
    exit 1
fi

# Calculate expected number of jobs
NUM_SAMPLES=${#SAMPLES[@]}
INTERVAL_DIR="/ceph/sharedfs/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/${PROJECT_NAME}/calling/interval/${REF_NAME}"
SCATTER_COUNT=$(find "${INTERVAL_DIR}" -maxdepth 1 -type f -name '*-scattered.interval_list' | wc -l)
EXPECTED_JOBS=$((NUM_SAMPLES * SCATTER_COUNT))

# Get the number of jobs from SLURM environment variable
ACTUAL_JOBS=${SLURM_ARRAY_TASK_COUNT:-$EXPECTED_JOBS} # Default for local testing

if [ "$EXPECTED_JOBS" -ne "$ACTUAL_JOBS" ]; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    echo "!!! ERROR: MISMATCH IN JOB ARRAY SIZE                      !!!" >&2
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    echo "Hardcoded array size in this script is: ${ACTUAL_JOBS} jobs." >&2
    echo "Calculated expected size is: ${EXPECTED_JOBS} jobs." >&2
    echo "  - Samples for ${REF_NAME}: ${NUM_SAMPLES}" >&2
    echo "  - Intervals: ${SCATTER_COUNT}" >&2
    echo "Please update the --array directive in this file ($0) and resubmit." >&2
    exit 1
fi
# Calculate which sample and interval this task should process
SAMPLE_IDX=$(( SLURM_ARRAY_TASK_ID / SCATTER_COUNT ))
INTERVAL_IDX=$(( SLURM_ARRAY_TASK_ID % SCATTER_COUNT ))

CRAM_FILE_TO_PROCESS="${SAMPLES[$SAMPLE_IDX]}"
INTERVAL_NUM=$(printf "%04d" "${INTERVAL_IDX}")
INTERVAL_FILE="${INTERVAL_DIR}/${INTERVAL_NUM}-scattered.interval_list"

echo "Executing generic HaplotypeCaller script for:"
echo "  - CRAM File: ${CRAM_FILE_TO_PROCESS}"
echo "  - Interval Index: ${INTERVAL_IDX}"

bash "${GENERIC_SCRIPT_PATH}" "${OUTPUT_BASE_DIR}" "${CRAM_FILE_TO_PROCESS}" "${REF_GENOME}" "${INTERVAL_FILE}" "${CORE}" "${PLOIDY}"

echo "Job task ${SLURM_ARRAY_TASK_ID} finished."
