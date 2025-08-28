#!/usr/bin/bash
#SBATCH --job-name=genome_size          # Job name for SLURM
#SBATCH --output=logs/%x_%j.log # Standard output log file
#SBATCH --error=logs/%x_%j.log  # Standard error log file
#SBATCH --cpus-per-task=1               # Number of CPU cores per task
#SBATCH --mem=4G                        # Memory per task
#SBATCH --time=00:30:00                 # Time limit for the job (HH:MM:SS) - Adjust as needed
#SBATCH --partition=intel-g4-al9_short_serial # SLURM partition to use

# written by BoHan Hou

# Set strict error handling
set -euo pipefail

# Create log directory if it doesn't exist
mkdir -p logs

echo "Starting genome size calculation at $(date)"
echo "----------------------------------------"

# Load modules (only once)
echo "Loading modules..."
module use /ceph/work/abrchmc/Software/modules
module load ucsc_tools/250521
echo "Modules loaded."
echo "----------------------------------------"

# Define an array of your reference FASTA file paths
declare -a REFERENCES=(
    "/ceph/work/abrchmc/Reference/mac/Musa_acuminata/DH-Pahang_v4.3/Musa_acuminata_pahang_v4.genome.fasta.gz"
	"/ceph/work/abrchmc/Reference/mac/Musa_acuminata/Baxijiao-haplotypes_Huang2023/Baxijiao.assembly.fasta.gz"
    "/ceph/work/abrchmc/Reference/mac/Musa_acuminata/Baxijiao-haplotypes_Li2023/Cavendish_chromosome.fasta.gz"
)

# Loop through each reference file in the array
for Ref in "${REFERENCES[@]}"; do
    Output=$(dirname "$Ref")

    echo "Processing reference: $Ref"
    echo "Output directory for this reference: $Output"

    # Change to the output directory for the current reference
    mkdir -p "$Output" # Ensure output directory exists
    cd "$Output" || { echo "ERROR: Could not change directory to $Output. Exiting."; exit 1; }
    echo "Changed directory to $Output"

    echo "Computing Effective Genome Size for $(basename "$Ref")..."
    # faCount: Counts bases, N's, and gaps in a FASTA file
    # awk -F'\t' '{print $1"\t"$2-$7}':
    #   -F'\t': Sets tab as field separator.
    #   $1: chromosome name
    #   $2: size of the sequence
    #   $7: number of N's (unknown bases)
    #   $2-$7: Calculates sequence length minus N's to get effective size
    faCount "$Ref" | awk -F'\t' '{print $1"\t"$2-$7}' > effectiveGenomeSize.txt
    echo "Effective genome size for $(basename "$Ref") computed and saved to effectiveGenomeSize.txt"
    echo "----------------------------------------"
done

echo "All effective genome sizes computed. Job finished at $(date)"
