#!/usr/bin/bash
# written by Puyam Tondonba Singh

input=/home/tony/project/GCTCV/trimmomatic/TRIM/ex
files=($input/*_R1_trimmed.fastq.gz)  # Use array to handle filenames with spaces
STAR_REF=/home/public_tools/NAS/share/Reference/Musa_acuminata/STAR_musa_v4_CorrectedMaACA8
output=/home/public_tools/NAS/share/tony/GCTCV/STAR/Align_sorted
output_local="/tmp"  # Local directory for temporary files
mkdir -p $output
cd $output

for file in "${files[@]}"
do
  echo "Processing file: $file"
  ID=$(basename "$file" | sed 's/_R1_trimmed.fastq.gz//')
  echo "Sample ID: $ID"

  # Check if both files exist before running STAR
  if [[ ! -f "${input}/${ID}_R1_trimmed.fastq.gz" || ! -f "${input}/${ID}_R2_trimmed.fastq.gz" ]]; then
    echo "ERROR: One or both read files are missing for $ID"
    continue
  fi

  echo "Start to align $ID..."
  echo "Using reads: ${input}/${ID}_R1_trimmed.fastq.gz and ${input}/${ID}_R2_trimmed.fastq.gz"

  # Run STAR with local temp directory for better file handling
  STAR \
    --runMode alignReads \
    --limitBAMsortRAM 9500000000 \
    --runThreadN 20 \
    --genomeDir $STAR_REF \
    --readFilesIn "${input}/${ID}_R1_trimmed.fastq.gz" "${input}/${ID}_R2_trimmed.fastq.gz" \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 1 \
    --alignIntronMax 1000 \
    --alignEndsType Local \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattrRGline ID:$ID LB:$ID PL:Illumina PU:$ID SM:$ID \
    --twopassMode Basic \
    --outFileNamePrefix ${output_local}/${ID}_ \
    --outTmpDir /tmp/${ID}_STARtmp

  # Check if STAR alignment was successful
  if [[ -f "${output_local}/${ID}_Aligned.sortedByCoord.out.bam" ]]; then
    # Move the output BAM file to the NAS
    mv ${output_local}/${ID}_Aligned.sortedByCoord.out.bam $output/
    samtools index -@ 30 ${output}/${ID}_Aligned.sortedByCoord.out.bam

    # Clean up temporary files and logs
    rm -f "$output_local/${ID}_Log.out" "$output_local/${ID}_Log.progress.out" "$output_local/${ID}_SJ.out.tab"
    echo "Alignment done for $ID."
  else
    echo "ERROR: STAR failed for $ID."
  fi
done

echo "All alignments completed."


/home/public_tools/NAS/share/Reference/Musa_acuminata/
