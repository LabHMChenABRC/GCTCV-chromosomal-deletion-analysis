#!/usr/bin/bash
#SBATCH --job-name=DHbamCoverage       # Job name for SLURM
#SBATCH --output=logs/DH/coverage_%A_%a.log # Standard output log file (%A for Job ID, %a for Array Task ID)
#SBATCH --error=logs/DH/coverage_%A_%a.log  # Standard error log file
#SBATCH --array=0-5                  # !! IMPORTANT: Adjust this based on the number of samples in Huang.sample.list (N-1) !!
#SBATCH --cpus-per-task=20           # Number of CPU cores per task
#SBATCH --mem=12G                    # Memory per task
#SBATCH --time=02:00:00              # Time limit for the job (HH:MM:SS) - Adjust as needed
#SBATCH --partition=edr1-al9_short   # SLURM partition to use

# written by BoHan Hou

# --- Set strict error handling for the submission script ---
set -euo pipefail

# --- User-defined variables for the array job ---
REF="/ceph/work/abrchmc/Reference/mac/Musa_acuminata/DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta.gz"
SAMPLE_LIST="${PWD}/DH.sample.list" # Path to your sample list file
OUTPUT_ROOT_DIR="/ceph/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/DH/coverage"

# DeepTools parameters
BIN_SIZE=10000
MAPQ_VALUES=(0) 
FLAG_EXCLUDE=2304 # not primary alignment(256); supplementary alignment(2048)

# --- Create logs directory if it doesn't exist ---
mkdir -p logs

echo "--- Starting coverage array job submission ---"

# Get the specific CRAM file path for this array task
CRAM_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")

# Extract ID (basename without .cram extension)
ID=$(basename "$CRAM_FILE" .cram)

# Determine the directory where effectiveGenomeSize.txt is located.
# This assumes effectiveGenomeSize.txt is in the same directory as the REF fasta.
EFFECTIVE_SIZE_DIR=$(dirname "$REF")
EFFECTIVE_SIZE_FILE="${EFFECTIVE_SIZE_DIR}/effectiveGenomeSize.txt"

# Core count from SLURM environment
CORE_COUNT="${SLURM_CPUS_PER_TASK:-1}" # Fallback to 1 if not set by Slurm

# --- Call the per-sample script ---
echo "Calling coverage.one.sh for sample: ${ID} (Task ID: ${SLURM_ARRAY_TASK_ID})"
echo "CRAM_FILE: ${CRAM_FILE}"
echo "EFFECTIVE_SIZE_FILE: ${EFFECTIVE_SIZE_FILE}"
echo "CORE_COUNT: ${CORE_COUNT}"
echo "BIN_SIZE: ${BIN_SIZE}"
echo "MAPQ_VALUES: ${MAPQ_VALUES[@]}"
echo "FLAG_EXCLUDE: ${FLAG_EXCLUDE}"
echo "OUTPUT_ROOT_DIR: ${OUTPUT_ROOT_DIR}"
echo "Sample ID: ${ID}"

# Execute the per-sample script with all necessary parameters
# Use 'bash' to explicitly run it, ensuring compatibility
for MAPQ_I in ${MAPQ_VALUES[@]}
do
  echo "  Processing with MAPQ: ${MAPQ_I}"
  bash coverage.one.sh \
    "$CRAM_FILE" \
    "$EFFECTIVE_SIZE_FILE" \
    "$CORE_COUNT" \
    "$BIN_SIZE" \
    "$MAPQ_I" \
    "$FLAG_EXCLUDE" \
    "$ID" \
    "$OUTPUT_ROOT_DIR"
done

echo "--- Finished coverage array job for task ${SLURM_ARRAY_TASK_ID} ---"
