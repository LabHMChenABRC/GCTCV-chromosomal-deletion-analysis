#!/usr/bin/bash
# This script is called by coverage.sbatch for each array task.
# written by BoHan Hou

# Set strict error handling
set -euo pipefail

# --- Receive arguments from the sbatch script ---
CRAM=$1                 # Full path to the input CRAM file
EFFECTIVE_SIZE_FILE=$2  # Full path to effectiveGenomeSize.txt
CORE=$3                 # Number of CPU cores to use for bamCoverage
BIN_SIZE=$4             # Bin size for bamCoverage
MAPQ=$5                 # Minimum mapping quality
FLAG_EXCLUDE=$6         # SAM flags to exclude
ID=$7                   # Sample ID (basename of CRAM file)
OUTPUT_ROOT_DIR=$8      # NEW: Root directory for BigWig output

echo "--- coverage.one.sh started for ${ID} ---"
echo "Received parameters:"
echo "  CRAM: ${CRAM}"
echo "  Effective Size File: ${EFFECTIVE_SIZE_FILE}"
echo "  Cores: ${CORE}"
echo "  Bin Size: ${BIN_SIZE}"
echo "  Min MapQ: ${MAPQ}"
echo "  Flag Exclude: ${FLAG_EXCLUDE}"
echo "  Sample ID: ${ID}"
echo "  Output Root Dir: ${OUTPUT_ROOT_DIR}" # Echo the new output directory
echo "----------------------------------------"

# --- Load modules ---
# Modules are loaded here as each array task is a fresh environment
echo "Loading modules..."
module use /ceph/work/abrchmc/Software/modules
module load deeptools
echo "Modules loaded."
echo "----------------------------------------"

# --- Validate effectiveGenomeSize.txt and get value ---
echo "Checking for effectiveGenomeSize.txt at: ${EFFECTIVE_SIZE_FILE}"
if [ ! -f "$EFFECTIVE_SIZE_FILE" ]; then
    echo "ERROR: effectiveGenomeSize.txt not found at $EFFECTIVE_SIZE_FILE. Please ensure it's generated."
    exit 1
fi

EFFECTIVE_GENOME_SIZE=$(grep total "$EFFECTIVE_SIZE_FILE" | awk '{print $2}')
if [ -z "$EFFECTIVE_GENOME_SIZE" ]; then
    echo "ERROR: 'total' entry not found or empty in effectiveGenomeSize.txt."
    exit 1
fi
echo "Using effective genome size: ${EFFECTIVE_GENOME_SIZE}"
echo "----------------------------------------"

# --- Main processing step: Run bamCoverage ---
echo "$(date) ... Start to calculate coverage with -p ${CORE}"

if [ -f "$CRAM" ]; then
    # Construct the full output path for the BigWig file
    # You might want to create subdirectories within OUTPUT_ROOT_DIR based on CRAM_FILE's original path structure
    # For now, let's just put it directly in OUTPUT_ROOT_DIR
    OUTPUT_BIGWIG="${OUTPUT_ROOT_DIR}/${ID}.RPGC.bin${BIN_SIZE}.MapQ${MAPQ}.bw"

    # Ensure the output directory exists before writing
    mkdir -p "$(dirname "$OUTPUT_BIGWIG")"

    bamCoverage \
        -b "$CRAM" \
        -o "$OUTPUT_BIGWIG" \
        -p "$CORE" \
        --exactScaling \
        --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" \
        --minMappingQuality "$MAPQ" \
        --normalizeUsing RPGC \
        --samFlagExclude "$FLAG_EXCLUDE" \
        --binSize "$BIN_SIZE"

    echo "bamCoverage completed for ${ID}. Output: ${OUTPUT_BIGWIG}"
else
    echo "ERROR: CRAM file ${CRAM} was not found, skipping it."
    exit 1 # Exit the array task if the input file is missing
fi

echo "--- coverage.one.sh finished for ${ID} ---"
