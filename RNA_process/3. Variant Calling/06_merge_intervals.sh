#!/usr/bin/bash
# written by Puyam Tondonba Singh

# Define directories and reference
input_dir="/home/tony/project/pnas_test/vcf/chr05/"
output_dir="/home/tony/project/pnas_test/vcf/chr05/merge"
mkdir -p "$output_dir"
fasta="/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta"

# Loop over unique sample files (using the .0000.g.vcf.gz file to extract sample name)
for file in "$input_dir"/*.marked.dup.cigarN.chr05.0000.g.vcf.gz; do
    # Extract sample name by removing the suffix
    sample=$(basename "$file" .marked.dup.cigarN.chr05.0000.g.vcf.gz)
    echo "Processing sample: $sample"
    
    # Create an input list for this sample
    sample_input_list="${output_dir}/${sample}_input.list"
    ls "${input_dir}/${sample}.marked.dup.cigarN.chr05."*.g.vcf.gz | sort > "$sample_input_list"
    echo "Input list for $sample created at: $sample_input_list"
    
    # Build the -V arguments from the input list
    vcf_args=""
    while IFS= read -r line; do
        vcf_args="$vcf_args -V $line"
    done < "$sample_input_list"
    
    # Define the merged output file for the sample
    merged_gvcf="${output_dir}/${sample}.merged.g.vcf.gz"
    
    echo "Merging files for sample $sample..."
    
    gatk --java-options "-Xmx16G" CombineGVCFs \
         -R "$fasta" \
         $vcf_args \
         -O "$merged_gvcf" \
         > "${output_dir}/${sample}_combine.log" 2>&1
    
    echo "Finished merging sample: $sample. Merged gVCF is available at: $merged_gvcf"
done

echo "All samples have been merged."
