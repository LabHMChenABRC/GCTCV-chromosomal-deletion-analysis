#!/bin/bash

# --- SLURM Configuration for Per-Chromosome Hard Filtering - Project: DHv4, All chromosomes ---
# This script runs the HF.one.sh pipeline for each chromosome in parallel.
# written by BoHan Hou
#
# Total Jobs = 1 REF * 11 chromosomes = 11. Array is 0-10.
# NOTE: The --array upper bound must match the number of chromosomes in the reference minus one.
#
# Before submitting, ensure the log directory exists:
# mkdir -p log/HF/DH_per_chr
#
#SBATCH --job-name=HF.DH.per_chr
#SBATCH --output=log/HF/DH_per_chr/%x_%A_%a.log
#SBATCH --error=log/HF/DH_per_chr/%x_%A_%a.log
#SBATCH --array=0-10 # Corresponds to 11 chromosomes
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=edr1-al9_short_serial

set -euo pipefail

# --- Configuration ---
PROJECT_NAME="DH"
REF_NAMES=("DHv4") # Add other reference names here, e.g., ("DHv4" "DHv2")
CORE_SCRIPT_PATH="./HF.one.sh"  # Path to the new per-chromosome filtering script
NUM_CHRS_PER_REF=11 # This must match the number of chromosomes in each reference
RAM_MB="${SLURM_MEM_PER_NODE}M" # SLURM_MEM_PER_NODE for memory in MB, e.g., 4096 for --mem=4G
CORE=${SLURM_CPUS_PER_TASK}

# --- Reference Genome Paths (Associative Array) ---
# This allows the script to dynamically select the correct reference path.
declare -A REF_GENOME_PATHS
REF_GENOME_PATHS["DHv4"]="/ceph/work/abrchmc/Reference/mac/Musa_acuminata/DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta.gz"
# Add paths for other references here, matching the names in REF_NAMES.
# REF_GENOME_PATHS["DHv2"]="/path/to/your/DHv2.fasta.gz"


# --- Path Definitions ---
PROJECT_DIR="/ceph/sharedfs/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/${PROJECT_NAME}"

REF_IDX=$(( SLURM_ARRAY_TASK_ID / NUM_CHRS_PER_REF ))
REF_NAME=${REF_NAMES[$REF_IDX]}

REF_GENOME=${REF_GENOME_PATHS[$REF_NAME]}
if [ -z "$REF_GENOME" ]; then
    echo "Error: No reference genome path defined for REF_NAME='${REF_NAME}' in REF_GENOME_PATHS." >&2
    exit 1
fi
REF_FAI="${REF_GENOME}.fai"

INPUT_GT_DIR="${PROJECT_DIR}/calling/GT/${REF_NAME}"
OUTPUT_HF_DIR="${PROJECT_DIR}/calling/HF/${REF_NAME}"

# --- Get Chromosome Name for this Task ---
# This section dynamically determines the chromosome for the current SLURM array task.
if [ ! -f "$REF_FAI" ]; then
    echo "Error: Reference index file not found: $REF_FAI" >&2
    exit 1
fi
mapfile -t CHR_NAMES < <(awk '{print $1}' "$REF_FAI")
CHR_IDX=$(( SLURM_ARRAY_TASK_ID % NUM_CHRS_PER_REF ))

if [ "$CHR_IDX" -ge "${#CHR_NAMES[@]}" ]; then
    echo "Error: SLURM array task ID ($SLURM_ARRAY_TASK_ID) resulted in an out-of-range chromosome index ($CHR_IDX)." >&2
    echo "This script is configured for ${#CHR_NAMES[@]} chromosomes per reference (indices 0-$((${#CHR_NAMES[@]} - 1)))." >&2
    echo "Please ensure the --array directive in the sbatch header matches the number of chromosomes (${#CHR_NAMES[@]})." >&2
    exit 1
fi
CHR_TO_PROCESS=${CHR_NAMES[$CHR_IDX]}

echo "--- Task Info ---"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Project: ${PROJECT_NAME}"
echo "REF Name: ${REF_NAME}"
echo "Chromosome: ${CHR_TO_PROCESS}"
echo "Input GT Dir: ${INPUT_GT_DIR}"
echo "Output HF Dir: ${OUTPUT_HF_DIR}"
echo "Cores: ${CORE}"
echo "RAM for GATK: ${RAM_MB}"
echo "-----------------"

# --- Execute the generic script ---
echo "Starting Per-Chromosome Hard Filtering task..."

bash "${CORE_SCRIPT_PATH}" \
    "${REF_GENOME}" \
    "${INPUT_GT_DIR}" \
    "${OUTPUT_HF_DIR}" \
    "${CHR_TO_PROCESS}" \
    "${CORE}" \
    "${RAM_MB}" || {
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    echo "!!! ERROR: Core script HF.one.sh failed.                   !!!" >&2
    echo "!!! Task: Hard Filtering for ${REF_NAME}/${CHR_TO_PROCESS} !!!" >&2
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >&2
    exit 1
}

echo "Job task ${SLURM_ARRAY_TASK_ID} for ${REF_NAME}/${CHR_TO_PROCESS} finished."
