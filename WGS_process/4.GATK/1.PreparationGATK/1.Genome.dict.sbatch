#!/usr/bin/bash
#SBATCH --job-name=mk.dict             # Job name for SLURM
#SBATCH --output=log/dict/%x_%A_%a.log # Standard output log file for the array job: JobName_ArrayJobID_ArrayTaskID.log
#SBATCH --error=log/dict/%x_%A_%a.log  # Standard error log file for the array job: JobName_ArrayJobID_ArrayTaskID.log
#SBATCH --array=0-6                    # Array for 7 reference files
#SBATCH --cpus-per-task=2              # CreateSequenceDictionary is single-threaded, 2 is plenty
#SBATCH --mem=4G                       # Memory per task
#SBATCH --time=0:30:00                 # 30 minutes is more than enough
#SBATCH --partition=intel-g4-al9_short # SLURM partition to use

# written by BoHan Hou

set -euo pipefail

# === Environment Setup ===
module use /ceph/work/abrchmc/Software/modules
module load gatk

# === Reference Genomes ===
# List of unique reference genomes to process
REFs=(
    "/ceph/work/abrchmc/Reference/mac/Musa_acuminata/DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta.gz"
)

# === Select Reference for this Array Task ===
# Check if the array task ID is valid
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#REFs[@]}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of bounds for REFs array size ${#REFs[@]}." >&2
    exit 1
fi
REF="${REFs[$SLURM_ARRAY_TASK_ID]}"

# === Define Output Path ===
# The dictionary file should be placed next to the reference FASTA.
# GATK convention is to name it by replacing .fasta.gz with .dict
OUTPUT_DICT="${REF%.fasta.gz}.dict"

# === Main Logic ===
echo "SLURM_JOB_ID: $SLURM_JOB_ID, Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing reference: $REF"
echo "Output dictionary: $OUTPUT_DICT"

if [[ -f "$OUTPUT_DICT" ]]; then
    echo "Dictionary already exists, skipping: $OUTPUT_DICT"
    exit 0
fi

echo "Creating dictionary for $REF..."
gatk --java-options "-Xmx4G" CreateSequenceDictionary \
    -R "$REF" \
    -O "$OUTPUT_DICT"

echo "Successfully created dictionary: $OUTPUT_DICT"
