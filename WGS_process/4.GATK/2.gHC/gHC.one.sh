#!/usr/bin/bash

# This is a generic GATK HaplotypeCaller script.

# --- Argument Parsing ---
if [ "$#" -ne 6 ]; then
    echo "Error: Missing arguments." >&2
    echo "Usage: $0 <OutputDir> <CramFilePath> <RefGenomeFastaGZip> <IntervalFile> <HMMthreads> <Ploidy>" >&2
    echo "Example: $0 /path/to/output /path/to/sample.cram /path/to/genome.fasta.gz /path/to/0000-scattered.interval_list 2 3" >&2
    exit 1
fi

OUTPUT_BASE_DIR=$1
CRAM_FILE=$2
REF_GENOME=$3
INTERVAL_FILE=$4
HMM_THREADS=$5
PLOIDY=$6


# --- Script Body ---
set -euo pipefail

echo "--- Job Start ---"
echo "CRAM File: ${CRAM_FILE}"

# --- Environment Setup ---
module use /ceph/work/abrchmc/Software/modules
module load gatk

# --- Task-specific File Selection ---
SAMPLE_ID=$(basename "$CRAM_FILE" .cram)
INTERVAL_NUM=$(basename "$INTERVAL_FILE" -scattered.interval_list)

# Define output file path
OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_ID}"
OUTPUT_GVCF="${OUTPUT_DIR}/${SAMPLE_ID}.${INTERVAL_NUM}.g.vcf.gz"

echo "Reference: ${REF_GENOME}"
echo "Sample ID: ${SAMPLE_ID}"
echo "Interval File: ${INTERVAL_FILE}"
echo "Output gVCF: ${OUTPUT_GVCF}"
echo "-------------------"

# --- Pre-run Checks ---
if [ ! -f "$REF_GENOME" ]; then echo "Error: Reference genome not found: $REF_GENOME"; exit 1; fi
if [ ! -f "${REF_GENOME}.fai" ]; then echo "Error: Reference index (.fai) not found for $REF_GENOME"; exit 1; fi
if [ ! -f "${REF_GENOME%.fasta.gz}.dict" ]; then echo "Error: Reference dictionary (.dict) not found for $REF_GENOME"; exit 1; fi
if [ ! -f "$CRAM_FILE" ]; then echo "Error: CRAM file not found: $CRAM_FILE"; exit 1; fi
if [ ! -f "$INTERVAL_FILE" ]; then echo "Error: Interval file not found: $INTERVAL_FILE"; exit 1; fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Idempotency check: skip if output and its index already exist
if [ -f "$OUTPUT_GVCF" ] && [ -f "${OUTPUT_GVCF}.tbi" ]; then
    echo "Output file and index already exist, skipping task."
    echo "--- Job End (Skipped) ---"
    exit 0
fi

# --- Main GATK Command ---
echo "Running GATK HaplotypeCaller..."
gatk --java-options "-Xmx14G -Djava.io.tmpdir=$OUTPUT_DIR/" HaplotypeCaller \
    -R "$REF_GENOME" \
    -I "$CRAM_FILE" \
    -L "$INTERVAL_FILE" \
    -O "$OUTPUT_GVCF" \
    --native-pair-hmm-threads "${HMM_THREADS}" \
    --sample-ploidy "${PLOIDY}" \
    -ERC BP_RESOLUTION \
    --linked-de-bruijn-graph \
    --output-mode EMIT_ALL_ACTIVE_SITES \
    --max-reads-per-alignment-start 0 \
    --dont-use-soft-clipped-bases true \
    --seconds-between-progress-updates 60

echo "HaplotypeCaller finished successfully."

echo "Indexing the output gVCF..."
gatk IndexFeatureFile -I "$OUTPUT_GVCF"

echo "Indexing finished successfully."
echo "--- Job End ---"
