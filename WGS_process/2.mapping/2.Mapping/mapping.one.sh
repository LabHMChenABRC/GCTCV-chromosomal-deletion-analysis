#!/usr/bin/env bash
# Exit immediately if a command exits with a non-zero status.
# Treat unset variables as an error when substituting.
# The exit status of a pipeline is the exit status of the last command that exited with a non-zero status.
set -euo pipefail

# written by BoHan Hou

# === Input Parameters ===
# Sample ID (e.g., 119, FM, TC5-2)
ID="$1"
# Path to the file containing sample information (R1, R2, RG string)
FILE_LIST="$2"
# Path to the reference genome
REF="$3"
# Number of CPU cores to use for BWA-MEM2 and Samtools
CORE="$4"
# Directory where output CRAM files will be stored
OUTPUT_DIR="$5"

echo "Running mapping for sample: $ID"
echo "Using file list: $FILE_LIST"
echo "Reference: $REF"
echo "Cores: $CORE"
echo "Output directory: $OUTPUT_DIR"

# === Extract FASTQ paths and Read Group (RG) string from FILE_LIST ===
# Grep for the line matching the current sample ID.
# Cut -f2 for R1 path, -f3 for R2 path, -f4- for the rest of the line (RG string).
# Using 'head -n 1' to ensure only the first match is taken if IDs are not strictly unique.
# This assumes the FILE_LIST format: ID <tab> R1_PATH <tab> R2_PATH <tab> @RG ...
R1=$(grep -P "^${ID}\t" "$FILE_LIST" | cut -f2 | head -n 1)
R2=$(grep -P "^${ID}\t" "$FILE_LIST" | cut -f3 | head -n 1)
RG=$(grep -P "^${ID}\t" "$FILE_LIST" | cut -f4 | head -n 1)

# Check if R1, R2, and RG were successfully extracted
if [ -z "$R1" ] || [ -z "$R2" ] || [ -z "$RG" ]; then
    echo "Error: Could not extract R1, R2, or RG information for sample ID: $ID from $FILE_LIST" >&2
    exit 1
fi

echo "R1: $R1"
echo "R2: $R2"
echo "RG: $RG"

# === Create Output Directory ===
# -p option prevents error if directory already exists
mkdir -p "$OUTPUT_DIR"

# === Change to Output Directory ===
# This simplifies output file paths, but full paths could also be used.
cd "$OUTPUT_DIR" || { echo "Error: Could not change to output directory $OUTPUT_DIR" >&2; exit 1; }

# === Perform Mapping and Sorting ===
# bwa-mem2 mem: Aligns reads to the reference genome.
#   -t $CORE: Specifies the number of threads.
#   -R "$RG": Adds the Read Group string to all reads. IMPORTANT for GATK.
#   $REF: Path to the reference genome.
#   $R1 $R2: Paths to the input FASTQ files (paired-end).
# | : Pipes the output (SAM format) of bwa-mem2 to samtools sort.
# samtools sort: Sorts the alignments by coordinate.
#   -@ $CORE: Specifies the number of threads for sorting.
#   -O cram: Output format is CRAM (compressed BAM, requires reference).
#   --reference $REF: Required for CRAM output.
#   -o $ID.cram: Output CRAM file name.
#   -: Reads input from stdin (from bwa-mem2).
echo "Starting BWA-MEM2 mapping and Samtools sorting..."
bwa-mem2 mem -t "$CORE" -R "$RG" "$REF" "$R1" "$R2" | \
samtools sort -@ "$CORE" -O cram --reference "$REF" -o "${ID}.cram" -

# Check if the mapping and sorting step was successful
if [ $? -ne 0 ]; then
    echo "Error: BWA-MEM2 mapping or Samtools sorting failed for sample $ID." >&2
    exit 1
fi
echo "Mapping and sorting completed: ${ID}.cram"

# === Index CRAM file and generate Flagstat report ===
# Check if the CRAM file was successfully created before indexing and flagstat.
if [ -f "${ID}.cram" ]; then
    echo "Indexing CRAM file..."
    samtools index -@ "$CORE" "${ID}.cram"
    if [ $? -ne 0 ]; then
        echo "Error: Samtools indexing failed for ${ID}.cram." >&2
        exit 1
    fi
    echo "CRAM index created: ${ID}.cram.crai"

    echo "Generating flagstat report..."
    samtools flagstat -@ "$CORE" "${ID}.cram" > "${ID}.cram.stat"
    if [ $? -ne 0 ]; then
        echo "Error: Samtools flagstat failed for ${ID}.cram." >&2
        exit 1
    fi
    echo "Flagstat report created: ${ID}.cram.stat"
    cat "${ID}.cram.stat"
else
    echo "Error: CRAM file ${ID}.cram not found after mapping and sorting. Skipping indexing and flagstat." >&2
    exit 1
fi

echo "Processing for sample $ID finished."
