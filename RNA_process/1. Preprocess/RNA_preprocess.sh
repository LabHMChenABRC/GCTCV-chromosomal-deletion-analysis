#!/usr/bin/bash
# written by Puyam Tondonba Singh

adapter=/home/public_tools/Software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa 
input=/home/tony/project/GCTCV/trimmomatic/trim_table
output=/home/tony/project/GCTCV/trimmomatic/TRIM/ex
logout=/home/tony/project/GCTCV/trimmomatic/TRIM/log_file
mkdir -p $output 
 
while IFS=$" " 
read -r ID SizeMin FileName1 FileName2 Path
do
    echo -n "$ID ... "
    trimmomatic-0.36 PE -threads 20 $Path$FileName1 $Path$FileName2 \
    $output/$ID"_R1_trimmed.fastq.gz" /dev/null $output/$ID"_R2_trimmed.fastq.gz" \
    /dev/null ILLUMINACLIP:$adapter:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:100 2>>${ID}.log
    echo "done"

done <$input/GCTCV.tsv 
echo "All finish"

