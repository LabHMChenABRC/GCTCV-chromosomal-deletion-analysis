#!/usr/bin/bash
#SBATCH --job-name=WGS_map_Huang        # Job name for SLURM
#SBATCH --output=log/Huang/%x_%A_%a.log          # Standard output log file
#SBATCH --error=log/Huang/%x_%A_%a.log           # Standard error log file
#SBATCH --array=0-5%3                   # Array job setup: 0-5 for 6 samples, %3 for max 3 concurrent jobs
#SBATCH --cpus-per-task=24              # Number of CPU cores per task
#SBATCH --mem=16G                       # Memory per task
#SBATCH --time=4:00:00                  # Time limit for the job (HH:MM:SS)
#SBATCH --partition=intel-g4-al9_short_serial # SLURM partition to use

# written by BoHan Hou

# === Environment Setup ===
# Load necessary modules for BWA-MEM2 and Samtools
module use /ceph/work/abrchmc/Software/modules
module load BWA-MEM
module load samtools

# === Script Options ===
# Path to the file containing sample information (IDs, R1, R2, RG string)
FILE_LIST="$PWD/WGS.RG.tsv"

# Output directory for mapped BAM/CRAM files.
# Creates a 'mapped_bams' directory relative to the current working directory.
#OUTPUT_DIR="/ceph/sharedfs/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/Baxijiao-haplotypes_Li2023/mapped_crams"
OUTPUT_DIR="/ceph/sharedfs/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/Baxijiao-haplotypes_Huang2023/mapped_crams"
LOG_DIR="$OUTPUT_DIR/log"
mkdir -p $LOG_DIR

# Reference genome path
#REF="/ceph/work/abrchmc/Reference/mac/Musa_acuminata/Baxijiao-haplotypes_Li2023/Cavendish_chromosome.fasta.gz"
REF="/ceph/work/abrchmc/Reference/mac/Musa_acuminata/Baxijiao-haplotypes_Huang2023/Baxijiao.assembly.fasta.gz"

# Number of cores to use per mapping job, derived from SLURM allocation
CORE="$SLURM_CPUS_PER_TASK"

# === Extract Sample ID for current array task ===
# Read all sample IDs from the first column of the FILE_LIST into a Bash array.
# The 'mapfile' (or 'readarray') command is used to read lines from a file into an array.
# 'cut -f1' extracts the first column (Sample ID).
mapfile -t IDS < <(cut -f1 "$FILE_LIST")

# Determine the sample ID for the current array task.
# SLURM_ARRAY_TASK_ID directly corresponds to the index in the array.
SAMPLE_ID="${IDS[$SLURM_ARRAY_TASK_ID]}"

# Check if the SAMPLE_ID was successfully retrieved
if [ -z "$SAMPLE_ID" ]; then
    echo "Error: Could not retrieve sample ID for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID." >&2
    exit 1
fi

echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Processing Sample ID: $SAMPLE_ID"
echo "Reference Genome: $REF"
echo "Output Directory: $OUTPUT_DIR"
echo "Cores per task: $CORE"
echo "Mapping Log file: $LOG_DIR/${SAMPLE_ID}.log"

# === Execute Mapping Script for the current sample ===
# Call the mapping.one.bash script with the required parameters.
# The parameters are: Sample ID, File List, Reference, Cores, Output Directory.
bash mapping.one.sh "$SAMPLE_ID" "$FILE_LIST" "$REF" "$CORE" "$OUTPUT_DIR" > "$LOG_DIR/${SAMPLE_ID}.log" 2>&1

# Check the exit status of the mapping.one.bash script
# The exit status is captured from the redirected command.
if [ $? -eq 0 ]; then
    echo "Mapping and processing for sample $SAMPLE_ID completed successfully. See $LOG_DIR/${SAMPLE_ID}.log for details."
else
    echo "Error: Mapping and processing for sample $SAMPLE_ID failed. Check $LOG_DIR/${SAMPLE_ID}.log for errors." >&2
    exit 1
fi
