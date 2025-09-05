#!/usr/bin/bash
# written by Puyam Tondonba Singh

fasta="/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta"
output=/home/tony/project/pnas_test/Intervals

for chr in $(cut -f 1 /home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta.fai | grep chr)

do 
  mkdir -p $output/$chr
  gatk SplitIntervals \
    -R $fasta \
    -L $chr \
    --scatter-count 20 \
    -O $output/$chr >/dev/null

done