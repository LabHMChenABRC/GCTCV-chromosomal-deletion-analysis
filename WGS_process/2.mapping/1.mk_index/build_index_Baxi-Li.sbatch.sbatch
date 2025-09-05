#!/usr/bin/bash
#SBATCH --job-name=bwa_index_Li_chr
#SBATCH --output=log/%x_%j.log
#SBATCH --error=log/%x_%j.log
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --partition=intel-g4-al9_short_serial

# written by BoHan Hou

# === Environment ===
module use /ceph/work/abrchmc/Software/modules
module load BWA-MEM
module load samtools
module load htslib
source $HOME/user_PATH
# === Parameters ===
genome=/ceph/work/abrchmc/Reference/mac/Musa_acuminata/Baxijiao-haplotypes_Li2023/Cavendish_genomic.fasta.gz

# === Log ===
TIMER=$(timer.sh start)
echo "Building BWA-MEM2 index for: $genome"
date "+%Y-%m-%d %H:%M:%S Start"

# === Prepare ===
genome_dir=$(dirname "$genome")
cd "$genome_dir" || exit 1

# Build .fai if missing
[ ! -f "${genome}.fai" ] && samtools faidx "$genome"

# Filter sequence IDs (confirm pattern)
cut -f1 "${genome}.fai" | sort -h | grep -P "Ban|Ze|Dh" > chrs_IDs.list

# Extract chromosomes
chromosome_file="$genome_dir/Cavendish_chromosome.fasta.gz"
samtools faidx "$genome" -r chrs_IDs.list | bgzip -@ 4 > "$chromosome_file"
samtools faidx "$chromosome_file"

# Indexing
bwa-mem2 index "$chromosome_file"

# === Done ===
date "+%Y-%m-%d %H:%M:%S Complete"
timer.sh end $TIMER
