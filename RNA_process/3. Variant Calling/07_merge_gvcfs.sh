!/usr/bin/bash
# written by Puyam Tondonba Singh

# Define directories and reference genome
input_dir="/home/tony/project/pnas_test/vcf/chr05/merge/"
output_dir="${input_dir}/combined"
fasta="/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta"

# Create output directory
mkdir -p "$output_dir"

# Collect all GVCF files (*.g.vcf.gz) in the input directory
gvcf_files=("$input_dir"/*.g.vcf.gz)

# Build the -V arguments for CombineGVCFs
vcf_args=$(printf ' -V %s' "${gvcf_files[@]}")

# Define the output file for the combined GVCF
combined_gvcf="${output_dir}/combined.chr05.g.vcf.gz"

# Run CombineGVCFs
gatk --java-options "-Xmx16G" CombineGVCFs \
     -R "$fasta" \
     $vcf_args \
     -O "$combined_gvcf" \
     > "$output_dir/combine_gvcfs.log" 2>&1

# Print completion message
echo "GVCF files combined successfully. The merged GVCF is available at: $combined_gvcf"
