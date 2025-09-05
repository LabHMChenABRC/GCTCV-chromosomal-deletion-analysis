#!/usr/bin/bash
# written by Puyam Tondonba Singh

bam=/home/tony/project/pnas_test/Align_sorted/*.bam
output_picard=/home/tony/project/pnas_test/Align_sorted/MARK_Dup
output_local="/tmp"  # Local directory for temporary files
mkdir -p $output_picard
cd $output_picard

for file in $bam
do
  echo $file  
  ID=$(basename $file | sed 's/_Aligned.sortedByCoord.out.bam//g')
  echo $ID
java -jar /home/public_tools/Software/picard/picard.jar \
      MarkDuplicates \
      I=$file \
      O=$output_picard/${ID}".marked_duplicates.bam" \
      M=$ID.txt
   samtools index -@ 30 ${ID}.marked_duplicates.bam

done
