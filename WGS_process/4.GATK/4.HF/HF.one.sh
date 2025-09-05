#!/usr/bin/bash
# written by BoHan Hou

set -euo pipefail

# --- Usage and Argument Parsing ---
START_TIME=$(date +%s)
if [[ $# -ne 6 ]]; then
    echo "Usage: $0 <reference.fasta.gz> <input_gt_dir> <output_dir> <chromosome> <cores> <ram>" >&2
    echo "  - reference.fasta.gz: Path to the reference genome." >&2
    echo "  - input_gt_dir:       Directory containing scattered 'genotyped.*.vcf.gz' files." >&2
    echo "  - output_dir:         Directory to save final filtered VCFs." >&2
    echo "  - chromosome:         The name of the chromosome to process (e.g., chr01)." >&2
    echo "  - cores:              Number of CPU cores to use." >&2
    echo "  - ram_gb:             Memory in GB for GATK Java options." >&2
    echo "Example: $0 ref.fa.gz /path/to/GT/DHv4 /path/to/HF/DHv4 chr01 2 4096M" >&2
    exit 1
fi

# --- Environment ---
module use /ceph/work/abrchmc/Software/modules
module load gatk
module load bcftools

# --- Parameters ---
REF=$1
INPUT_GT_DIR=$2
OUTPUT_DIR=$3
CHR=$4
CORE=$5
RAM_MB=$6
JAVA_OPTS="-Xmx${RAM_MB}"

# --- Main Workflow ---
echo "Starting Hard Filtering for Chromosome: ${CHR}"
echo "Input GT Directory: ${INPUT_GT_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Error: Could not change to output directory $OUTPUT_DIR"; exit 1; }

# --- Define filenames for this chromosome's intermediate and final files ---
GT_VCF_LIST="gt_vcf_for_${CHR}.list"
MERGED_CHR_VCF="00_merged.${CHR}.vcf.gz"
SNP_VCF="01_${CHR}_snps.vcf.gz"
INDEL_VCF="01_${CHR}_indels.vcf.gz"
SNP_HF_VCF="HF_${CHR}_snps.vcf.gz"
INDEL_HF_VCF="HF_${CHR}_indels.vcf.gz"

# === STEP 1: Merge GT files for a single chromosome ===
echo "STEP 1: Finding and merging GT files for ${CHR}..."
find "${INPUT_GT_DIR}" -type f -name "genotyped.*.vcf.gz" | sort > "$GT_VCF_LIST"
if [ ! -s "$GT_VCF_LIST" ]; then
    echo "ERROR: No 'genotyped.*.vcf.gz' files found in ${INPUT_GT_DIR}. Exiting." >&2
    exit 1
fi

# Use bcftools to merge and subset to the target chromosome in one go.
bcftools concat --threads "${CORE}" -f "${GT_VCF_LIST}" -r "${CHR}" -a -Oz -o "${MERGED_CHR_VCF}"

# Check if the merged file contains any variants. If not, the chromosome has no variants to filter.
if ! bcftools view -H "${MERGED_CHR_VCF}" > /dev/null; then
    echo "No variants found for chromosome ${CHR} in the input files. Skipping."
    rm -f "${GT_VCF_LIST}" "${MERGED_CHR_VCF}"
    exit 0
fi

bcftools index --threads "${CORE}" --tbi "${MERGED_CHR_VCF}"
echo "Successfully created merged VCF for ${CHR}: ${MERGED_CHR_VCF}"

# === STEP 2: Split into INDELs and SNPs ===
echo "STEP 2: Splitting SNPs and INDELs for ${CHR}..."
gatk --java-options "${JAVA_OPTS}" SelectVariants -R "${REF}" -V "${MERGED_CHR_VCF}" -select-type SNP -O "${SNP_VCF}"
gatk --java-options "${JAVA_OPTS}" SelectVariants -R "${REF}" -V "${MERGED_CHR_VCF}" -select-type INDEL -O "${INDEL_VCF}"

# === STEP 3: Hard filter SNPs and INDELs ===
echo "STEP 3: Applying hard filters..."
gatk --java-options "${JAVA_OPTS}" VariantFiltration \
    -R "${REF}" -V "${SNP_VCF}" -O "${SNP_HF_VCF}" \
    --verbosity ERROR \
    --cluster-window-size 25 --cluster-size 3 \
    --filter-name "QD2" -filter "QD < 2.0" \
    --filter-name "QUAL30" -filter "QUAL < 30.0" \
    --filter-name "SOR3" -filter "SOR > 3.0" \
    --filter-name "FS60" -filter "FS > 60.0" \
    --filter-name "MQ40" -filter "MQ < 40.0" \
    --filter-name "MQRankSum-12.5" -filter "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum-8" -filter "ReadPosRankSum < -8.0"

gatk --java-options "${JAVA_OPTS}" VariantFiltration \
    -R "${REF}" -V "${INDEL_VCF}" -O "${INDEL_HF_VCF}" \
    --verbosity ERROR \
    --filter-name "QD2" -filter "QD < 2.0" \
    --filter-name "QUAL30" -filter "QUAL < 30.0" \
    --filter-name "FS200" -filter "FS > 200.0" \
    --filter-name "ReadPosRankSum-20" -filter "ReadPosRankSum < -20.0"

# === STEP 4: Index results ===
echo "STEP 4: Verifying indexes for final files..."
if [ ! -f "${SNP_HF_VCF}.tbi" ]; then
    echo "Warning: Index for filtered SNPs not found. Creating manually."
    gatk IndexFeatureFile -I "${SNP_HF_VCF}"
fi
if [ ! -f "${INDEL_HF_VCF}.tbi" ]; then
    echo "Warning: Index for filtered INDELs not found. Creating manually."
    gatk IndexFeatureFile -I "${INDEL_HF_VCF}"
fi

# === STEP 5: Remove intermediate files ===
echo "STEP 5: Cleaning up intermediate files..."
rm -f "${GT_VCF_LIST}" \
      "${MERGED_CHR_VCF}" "${MERGED_CHR_VCF}.tbi" \
      "${SNP_VCF}" "${SNP_VCF}.tbi" \
      "${INDEL_VCF}" "${INDEL_VCF}.tbi"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(( (DURATION % 3600) / 60 ))
SECONDS=$((DURATION % 60))

echo "Hard Filtering for ${CHR} Completed Successfully."
echo "--------------------------------------------------"
echo "Total execution time: $(printf '%02d:%02d:%02d' $HOURS $MINUTES $SECONDS) (H:M:S)"
echo "--------------------------------------------------"
