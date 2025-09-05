#!/usr/bin/bash
# written by BoHan Hou

set -euo pipefail

module use /ceph/work/abrchmc/Software/modules
module load BBTools

# ARGVs:
#   $1 = sample ID
#   $2 = File_list（ID\tR1\tR2）
#   $3 = OUTPUT_DIR
#   $4 = RAM（eg. "8G"）
ID="$1"
FILE_LIST="$2"
OUTPUT_DIR="$3"
RAM="$4"

# Run Java programe BBTools with Max $RAM
JAVA_RAM="-Xmx${RAM}"

# Get R1 / R2
R1=$(grep -P "^${ID}\t" "$FILE_LIST" | cut -f2)
R2=$(grep -P "^${ID}\t" "$FILE_LIST" | cut -f3)

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "$(date '+%F %T') ... Process $ID"

echo "$(date '+%F %T') ... Remove duplicated reads"
clumpify.sh \
    in="$R1" in2="$R2" \
    out="${ID}.Dup.1P.fastq.gz" out2="${ID}.Dup.2P.fastq.gz" \
    $JAVA_RAM \
    dedupe reorder

echo "$(date '+%F %T') ... Quality trimming (qtrim=rl trimq=20 minlength=100)"
bbduk.sh \
    in="${ID}.Dup.1P.fastq.gz" in2="${ID}.Dup.2P.fastq.gz" \
    out="${ID}.1P.fastq.gz" out2="${ID}.2P.fastq.gz" \
    $JAVA_RAM \
    qtrim=rl trimq=20 minlength=100

# Remove temp fles after trimming
if [[ -f "${ID}.1P.fastq.gz" && -f "${ID}.2P.fastq.gz" ]]; then
    rm -f "${ID}.Dup.1P.fastq.gz" "${ID}.Dup.2P.fastq.gz"
fi

echo "$(date '+%F %T') ... Finish $ID"
