#!/usr/bin/bash
# written by Puyam Tondonba Singh

bam=/home/tony/project/pnas_test/Align_sorted/MARK_Dup/*.bam
fasta="/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta"
SplitNCigar="/home/tony/project/pnas_test/Align_sorted/SplitNcigar"
mkdir -p $SplitNCigar


for file in $bam
do
  ID=$(basename "$file" | sed "s/_duplicates.bam//g")
  output_file="$SplitNCigar/${ID}.dup.cigarN.bam"

  if [[ ! -f "$output_file" ]]; then
    echo "Processing $file"

    # Run GATK SplitNCigarReads (run sequentially)
    gatk SplitNCigarReads \
      -R "$fasta" \
      -I "$file" \
      -O "$output_file" \
      > "$SplitNCigar/${ID}.log" 2>&1

    # Index the BAM file with samtools
    samtools index -@ 30 "$output_file"
  else
    echo "Skipping $file (already processed)"
  fi

  # Uncomment the following line if you truly want ONLY ONE FILE processed:
  # break
done




