#!/usr/bin/bash
# written by Puyam Tondonba Singh

# Exit if any command fails.
set -e

# Define file paths
GTF="/home/public_tools/NAS/share/Reference/Musa_acuminata/Musa_acuminata_DH-Pahang_v4.3/Musa_acuminata_pahang_v4.CorrectedMaACA8.gtf"
GENOME_FASTA="/home/public_tools/NAS/share/Reference/Musa_acuminata/Musa_acuminata_DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta"
RSEM_REF_DIR="/home/tony/project/DEG_test/RSEM_ref/corrected/"

# Create output directory
mkdir -p "${RSEM_REF_DIR}"
cd "${RSEM_REF_DIR}"

# Extract gene-to-transcript mapping
awk -F'\t' '
    match($0, /gene_id "([^"]+)"/){gene=substr($0, RSTART+9, RLENGTH-10)} 
    match($0, /transcript_id "([^"]+)"/){transcript=substr($0, RSTART+15, RLENGTH-16); print gene "\t" transcript}
' "${GTF}" | sort -k1 -k2 | uniq > Gene.txt

# Run RSEM reference preparation
rsem-prepare-reference \
    --gtf "${GTF}" \
    --transcript-to-gene-map Gene.txt \
    --star \
    -p 30 \
    "${GENOME_FASTA}" \
    Musav4

echo "RSEM reference completed."