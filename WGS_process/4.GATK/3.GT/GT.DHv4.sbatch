#!/bin/bash

# --- SLURM Configuration for Genotyping - Project: DHv4, All genomic intervals ---
# This script performs joint-genotyping for each intervals in parallel.
# written by BoHan Hou
#
# Total Jobs = 1 REF * 30 intervals = 30. Array is 0-29.
#
# Before submitting, ensure the log directory exists:
# mkdir -p log/GT/DH/all
#
#SBATCH --job-name=GT.DH.all
#SBATCH --output=log/GT/DH/all/%x_%A_%a.log
#SBATCH --error=log/GT/DH/all/%x_%A_%a.log
#SBATCH --array=0-29
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=intel-g4-al9_short_serial

set -euo pipefail

# --- Configuration ---
PROJECT_NAME="DH"
REF_NAMES=("DHv4") # Order is important for mapping
CORE_SCRIPT_PATH="./GT.sh" # Path to the core genotyping script
PLOIDY=3 # Default ploidy
NUM_INTERVALS_PER_REF=30 # This must match the scatter count
RAM_MB="${SLURM_MEM_PER_NODE}M" # SLURM_MEM_PER_NODE for memory in MB, e.g., 4096 for --mem=4G
CORE=${SLURM_CPUS_PER_TASK}


# --- Map SLURM task ID to copy and interval ---
REF_IDX=$(( SLURM_ARRAY_TASK_ID / NUM_INTERVALS_PER_REF ))
INTERVAL_IDX=$(( SLURM_ARRAY_TASK_ID % NUM_INTERVALS_PER_REF ))
EX_INTERVAL_BED=$PWD/DHv4.XL.bed

REF_NAME=${REF_NAMES[$REF_IDX]}

# --- Path Definitions ---
PROJECT_DIR="/ceph/sharedfs/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/${PROJECT_NAME}"
REF_GENOME_PATHS["DHv4"]="/ceph/work/abrchmc/Reference/mac/Musa_acuminata/DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta.gz"
# Add paths for other references here, matching the names in REF_NAMES.
# REF_GENOME_PATHS["DHv2"]="/path/to/your/DHv2.fasta.gz"
REF_GENOME=${REF_GENOME_PATHS[$REF_NAME]}
if [ -z "$REF_GENOME" ]; then
    echo "Error: No reference genome path defined for REF_NAME='${REF_NAME}' in REF_GENOME_PATHS." >&2
    exit 1
fi

INTERVAL_DIR="${PROJECT_DIR}/calling/interval/${REF_NAME}"
GVCF_BASE_DIR="${PROJECT_DIR}/calling/gVCF/${REF_NAME}"
OUTPUT_DIR="${PROJECT_DIR}/calling/GT/${REF_NAME}"

echo "--- Task Info ---"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "REF Index: $REF_IDX -> ${REF_NAME}"
echo "Interval Index: $INTERVAL_IDX"
echo "-----------------"

# Get the list of all interval files for this haplotype
mapfile -t INTERVAL_FILES < <(find "${INTERVAL_DIR}" -maxdepth 1 -type f -name '*-scattered.interval_list' | sort)

# Sanity check the number of found interval files
if [ "${#INTERVAL_FILES[@]}" -ne "$NUM_INTERVALS_PER_REF" ]; then
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    echo "!!! ERROR: MISMATCH IN INTERVAL COUNT                                 !!!" >&2
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    echo "Expected ${NUM_INTERVALS_PER_REF} interval files, but found ${#INTERVAL_FILES[@]} in ${INTERVAL_DIR}" >&2
    echo "Please check the output of 2.mk.interval.sbatch or update NUM_INTERVALS_PER_REF in this script." >&2
    exit 1
fi

# Check that the calculated interval index is valid
if [ "$INTERVAL_IDX" -ge "${#INTERVAL_FILES[@]}" ]; then
    echo "Error: Calculated interval index (${INTERVAL_IDX}) is out of range for the number of files found (${#INTERVAL_FILES[@]})." >&2
    exit 1
fi

# Select the interval file for this specific task
INTERVAL_FILE_FOR_TASK="${INTERVAL_FILES[$INTERVAL_IDX]}"

# --- Execute the generic script ---
echo "Starting genotyping task for:"
echo "  - Project: ${PROJECT_NAME}"
echo "  - Reference: ${REF_NAME}"
echo "  - Interval File: ${INTERVAL_FILE_FOR_TASK}"

bash "${CORE_SCRIPT_PATH}" \
    "${REF_GENOME}" \
    "${INTERVAL_FILE_FOR_TASK}" \
    "${GVCF_BASE_DIR}" \
    "${OUTPUT_DIR}" \
    "${CORE}" \
    "${RAM_MB}" \
    "${PLOIDY}" \
    "${EX_INTERVAL_BED}"
    

echo "Job task ${SLURM_ARRAY_TASK_ID} finished."
