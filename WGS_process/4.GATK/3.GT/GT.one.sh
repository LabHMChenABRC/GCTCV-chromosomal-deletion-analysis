#!/usr/bin/bash
set -euo pipefail

if [[ $# -lt 6 || $# -gt 8 ]]; then
    echo "Usage: $0 <reference_file> <interval_file> <gVCF_base_dir> <output_dir> <num_cores> <RAM_in_GB> <Ploidy> [exclude_intervals_file]" >&2
    echo "Example: $0 ref.fa interval.list /path/to/gvcfs /path/to/output 4096M 8 3 [exclude.list]" >&2
    exit 1
fi

# --- Environment ---
module use /ceph/work/abrchmc/Software/modules
module load gatk

# --- Parameters ---
REF=$1
INT_FILE=$2
GVCF_BASE_DIR=$3
OUTPUT_DIR=$4
CORE=$5
RAM_MB=$6
PLOIDY=${7:-3} # Default ploidy is 3 if not provided
XL_FILE=${8:-} # Optional: 7th argument for exclusion intervals
JAVA_OPTS="-Xmx${RAM_MB} -XX:ParallelGCThreads=${CORE}"

IntervalNum="$(basename "$INT_FILE" -scattered.interval_list)"

# --- Output File Paths ---
mkdir -p "$OUTPUT_DIR"
GVCF_LIST_FILE="$OUTPUT_DIR/gvcfs_for_${IntervalNum}.list" # List of gVCFs for this interval
COMBINED_GVCF="$OUTPUT_DIR/combined.${IntervalNum}.g.vcf.gz"
FINAL_VCF="$OUTPUT_DIR/genotyped.${IntervalNum}.vcf.gz"

echo "--- Starting Genotyping for Interval: ${IntervalNum} ---"

# --- Step 1: Find all scattered gVCFs for this interval ---
echo "$(date "+%Y-%m-%d %H:%M:%S") ... Finding gVCFs in ${GVCF_BASE_DIR}"
find "${GVCF_BASE_DIR}" -type f -name "*\.${IntervalNum}\.g\.vcf\.gz" | sort > "$GVCF_LIST_FILE"

if [ ! -s "$GVCF_LIST_FILE" ]; then
    echo "Error: No gVCF files found for interval ${IntervalNum} in ${GVCF_BASE_DIR}" >&2
    exit 1
fi
echo "Found $(wc -l < "$GVCF_LIST_FILE") gVCF files for this interval."

# --- Step 2: Combine gVCFs ---
echo "$(date "+%Y-%m-%d %H:%M:%S") ... Combining gVCFs"
gatk --java-options "${JAVA_OPTS}" CombineGVCFs \
	-R "$REF" \
	-L "$INT_FILE" \
	--variant "$GVCF_LIST_FILE" \
	--seconds-between-progress-updates 60 \
    --output "$COMBINED_GVCF"

# --- Step 3: Joint-Genotype the combined gVCF ---
echo "$(date "+%Y-%m-%d %H:%M:%S") ... Genotyping combined gVCF"

# Build the GenotypeGVCFs command arguments in an array for flexibility
GATK_ARGS=(
	-R "$REF"
	-L "$INT_FILE"
	--sample-ploidy ${PLOIDY}
	--seconds-between-progress-updates 60
	-stand-call-conf 10
	--variant "$COMBINED_GVCF"
	--output "$FINAL_VCF"
)

# Conditionally add the exclude-intervals argument if the file path was provided
if [[ -n "$XL_FILE" ]]; then
    echo "Excluding intervals from: $XL_FILE"
    GATK_ARGS+=("--exclude-intervals" "$XL_FILE")
fi

gatk --java-options "${JAVA_OPTS}" GenotypeGVCFs "${GATK_ARGS[@]}"

# --- Step 4: Cleanup ---
echo "$(date "+%Y-%m-%d %H:%M:%S") ... Cleaning up intermediate files"
[ -f "$COMBINED_GVCF" ] && [ -f "${COMBINED_GVCF}.tbi" ] && [ -f "${FINAL_VCF}.tbi" ] && rm -f "$COMBINED_GVCF" "${COMBINED_GVCF}.tbi" "$GVCF_LIST_FILE"
echo "$(date "+%Y-%m-%d %H:%M:%S") ... Finished Interval: ${IntervalNum}"
