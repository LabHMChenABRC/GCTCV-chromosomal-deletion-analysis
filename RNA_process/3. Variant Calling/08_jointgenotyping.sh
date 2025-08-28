#!/usr/bin/bash
# written by Puyam Tondonba Singh

# Define directories and reference genome
input_dir="/home/tony/project/pnas_test/vcf/chr05/merge/combined"
output_dir="/home/tony/project/pnas_test/vcf/chr05/merge/joint"
fasta="/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta"

mkdir -p $output_dir

# Run GenotypeGVCFs
gatk GenotypeGVCFs \
  -R "$fasta" \
  -V "$input_dir/combined.chr05.g.vcf.gz" \
  -O "$output_dir/joint_genotyped.vcf.gz"
