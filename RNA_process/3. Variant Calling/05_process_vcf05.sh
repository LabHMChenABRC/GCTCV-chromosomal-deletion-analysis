#!/usr/bin/bash
# written by Puyam Tondonba Singh

# Define input and output directories
bam_dir="/home/public_tools/NAS/share/tony/GCTCV/STAR/SplitNCigar"
fasta="/home/tony/project/reference/Pahang_v4/Musa_acuminata_pahang_v4.fasta"
output_dir="/home/tony/project/GVCF/chr05/"
log_dir="/home/tony/project/GVCF/log"
interval_dir="/home/tony/project/GCTCV/Intervals/"

# Create output directories
mkdir -p $output_dir $log_dir

# Set max parallel jobs (adjust based on system resources)
MAX_JOBS=6  # Change this based on CPU cores available
job_count=0

# Process each BAM file
for bam in $bam_dir/*.bam; do
  sample_name=$(basename "$bam" .bam)
  
  for chr_id in chr05; do
    for interval in $interval_dir/$chr_id/*.interval_list; do
      interval_no=$(basename "$interval")
      interval_no=${interval_no%-scattered.interval_list}

      # Skip if output already exists
      if [ -f ${output_dir}/${sample_name}.${chr_id}.${interval_no}.g.vcf.gz ]; then
        echo "Skipping ${sample_name} ${chr_id} (already processed)"
        continue
      fi

      # Run HaplotypeCaller in GVCF mode (background process)
      gatk --java-options "-Xmx12G" HaplotypeCaller \
        -R $fasta \
        -I $bam \
        -L $interval \
        -ERC GVCF \
        -stand-call-conf 20.0 --max-reads-per-alignment-start 0 \
        -ploidy 3 \
        --dont-use-soft-clipped-bases \
        -O ${output_dir}/${sample_name}.${chr_id}.${interval_no}.g.vcf.gz \
        > ${log_dir}/${sample_name}.${chr_id}.${interval_no}.log 2>&1 &

      # Track the number of running jobs
      job_count=$((job_count + 1))

      # Limit jobs to prevent overloading the system
      if [[ $job_count -ge $MAX_JOBS ]]; then
        wait  # Wait for background jobs to finish before starting new ones
        job_count=0
      fi
    done
  done
done

# Final wait to ensure all background jobs are completed
wait

echo "GVCF generation completed!"
