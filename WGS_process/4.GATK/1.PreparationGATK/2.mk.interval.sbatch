#!/usr/bin/bash
#SBATCH --job-name=mk.interval
#SBATCH --output=log/interval/%x_%A_%a.log
#SBATCH --error=log/interval/%x_%A_%a.log
#SBATCH --array=0-5                    # Array for the 6 reference genomes
#SBATCH --cpus-per-task=2              # Low resource task
#SBATCH --mem=4G                       # Memory per task
#SBATCH --time=0:30:00                 # 30 minutes is sufficient
#SBATCH --partition=intel-g4-al9_short # SLURM partition to use

# written by BoHan Hou

set -euo pipefail

# === Environment Setup ===
module use /ceph/work/abrchmc/Software/modules
module load gatk

# === Configuration ===
# List of unique reference genomes to process
REFs=(
    "/ceph/work/abrchmc/Reference/mac/Musa_acuminata/DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta.gz"
)
# Corresponding short names for output directories
COPY_NAMES=("DHv4")

# Number of interval files to create per reference
SCATTER_COUNT=30

# === Select Task-Specific Files ===
# Check if the array task ID is valid
if [ -z "${SLURM_ARRAY_TASK_ID+x}" ]; then
    echo "Error: This script is designed to be run as a SLURM array job." >&2
    echo "Please set SLURM_ARRAY_TASK_ID for testing, e.g., export SLURM_ARRAY_TASK_ID=0" >&2
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" -ge "${#REFs[@]}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID ${SLURM_ARRAY_TASK_ID} is out of bounds for REFs array size ${#REFs[@]}." >&2
    exit 1
fi

REF="${REFs[$SLURM_ARRAY_TASK_ID]}"
COPY_NAME="${COPY_NAMES[$SLURM_ARRAY_TASK_ID]}"

# === Determine Output Directory ===
# Dynamically determine the base output directory from the reference path
if [[ "$REF" == *"pahang_v4"* ]]; then
    BASE_OUTPUT_DIR="/ceph/sharedfs/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/DH/calling/interval"
else 
    echo "Error: Cannot determine project (Huang2023/Li2023/pahang_v4) from reference path: $REF" >&2
    exit 1
fi
OUTPUT_DIR="$BASE_OUTPUT_DIR/$COPY_NAME/"

# === Main Logic ===
echo "--- Task ${SLURM_ARRAY_TASK_ID} ---"
echo "Reference: $REF"
echo "Copy Name: $COPY_NAME"
echo "Output Dir: $OUTPUT_DIR"

# Define path to the dictionary file, which should exist alongside the reference
REF_DICT="${REF%.fasta.gz}.dict"

if [[ ! -f "$REF_DICT" ]]; then
    echo "Error: Dictionary file not found: $REF_DICT" >&2
    echo "Please run 1.Genome.dict.sh first." >&2
    exit 1
fi

# Check if intervals already exist to avoid re-running
LAST_INTERVAL_FILE="$OUTPUT_DIR/$(printf "%04d" $((SCATTER_COUNT - 1)))-scattered.interval_list"
if [[ -f "$LAST_INTERVAL_FILE" ]]; then
    echo "Interval files already exist in $OUTPUT_DIR. Skipping."
    exit 0
fi

echo "Creating output directory: $OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Error: Failed to change directory to $OUTPUT_DIR"; exit 1; }

TEMP_INTERVAL_LIST="contigs_to_split.list"

echo "Extracting contig names from $REF_DICT..."
# Extract sequence names (SN) from the dictionary.
# The original command was incorrect. This awk command is more robust and correct.
# It filters for @SQ lines, then extracts and cleans the second field (e.g., SN:chr01.1 -> chr01.1).
awk -F'\t' '$1 == "@SQ" { sub("SN:","",$2); print $2 }' "$REF_DICT" > "$TEMP_INTERVAL_LIST"

if [ ! -s "$TEMP_INTERVAL_LIST" ]; then
    echo "Error: No contigs extracted from dictionary file. The temporary list is empty." >&2
    exit 1
fi

echo "Running GATK SplitIntervals to create $SCATTER_COUNT interval files..."
gatk SplitIntervals \
  -R "$REF" \
  -L "$TEMP_INTERVAL_LIST" \
  --scatter-count "$SCATTER_COUNT" \
  --subdivision-mode INTERVAL_SUBDIVISION \
  --output ./ \
  --interval-merging-rule OVERLAPPING_ONLY

rm "$TEMP_INTERVAL_LIST"

echo "Successfully created interval files in $OUTPUT_DIR"
echo "--- Task ${SLURM_ARRAY_TASK_ID} Finished ---"
