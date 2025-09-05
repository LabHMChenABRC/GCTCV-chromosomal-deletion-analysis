#!/usr/bin/bash
#SBATCH --job-name=bwa_index_Huang
#SBATCH --output=log/%x_%j.log
#SBATCH --error=log/%x_%j.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0:30:00
#SBATCH --partition=intel-g4-al9_short

# written by BoHan Hou

# === Environment ===
module use /ceph/work/abrchmc/Software/modules
module load BWA-MEM

# === Parameters ===
genome=/ceph/work/abrchmc/Reference/mac/Musa_acuminata/Baxijiao-haplotypes_Huang2023/Baxijiao.assembly.fasta.gz

# === Log ===
echo "Building BWA-MEM2 index for: $genome"
date "+%Y-%m-%d %H:%M:%S Start"

# === Indexing ===
bwa-mem2 index "$genome"

# === Done ===
date "+%Y-%m-%d %H:%M:%S Complete"
echo "Index files generated in: $index_dir"

