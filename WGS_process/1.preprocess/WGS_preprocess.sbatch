#!/usr/bin/bash
#SBATCH --job-name=WGS_preproc
#SBATCH --output=log/WGS_preproc_%x_%j.log
#SBATCH --error=log/WGS_preproc_%x_%j.log
#SBATCH --array=0-5%6           # 6 samples, run 6 jobs at once
#SBATCH --cpus-per-task=12
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --partition=intel-g4-al9_short_serial

# written by BoHan Hou

# === options ===
FILE_LIST=$PWD/sample.list
OUTPUT_DIR=/ceph/work/abrchmc/tcman/NGS/Banana_PNAS/WGS/preprocess
MAXJOB=6
RAM="${SLURM_MEM_PER_NODE}M"
SAMPLE_INDEX=$(( SLURM_ARRAY_TASK_ID % MAXJOB ))

# Get Sample IDs
mapfile -t IDS < <(cut -f1 "$FILE_LIST" )
ID=${IDS[$SAMPLE_INDEX]}
echo $ID
# remove dup reads and trim by quality score
bash WGS_preprocess.one.sh "$ID" "$FILE_LIST" "$OUTPUT_DIR" "$RAM"