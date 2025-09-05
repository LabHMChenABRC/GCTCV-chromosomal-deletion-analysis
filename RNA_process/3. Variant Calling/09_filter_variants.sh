#!/usr/bin/bash
# written by Puyam Tondonba Singh

# Define directories and reference genome
input_dir="/home/tony/project/pnas_test/vcf/chr05/merge/joint"
output_dir="/home/tony/project/pnas_test/vcf/chr05/merge/joint/filter"
fasta="/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta"

# Define input and output files
input_vcf="$input_dir/joint_genotyped.vcf.gz"
output_vcf="$output_dir/filtered_variants.vcf.gz"


mkdir -p $output_dir

# Run VariantFiltration to filter based on combined criteria
gatk VariantFiltration \
  -R "$fasta" \
  -V "$input_vcf" \
  --filter-name "PASS" \
  --filter-expression "QD < 2.0 ||FS > 60.0" \
  --cluster-size 3 \
  --cluster-window-size 35 \
  -O "$output_vcf"
