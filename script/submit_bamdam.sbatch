#!/bin/bash
#SBATCH --job-name=PhD26_bamdam
#SBATCH --output=/home/usr/shenmj/2026-PhD_project/tmp/logs/snakemake_bamdam_%j.log
#SBATCH --error=/home/usr/shenmj/2026-PhD_project/tmp/logs/snakemake_bamdam_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# ─── activate environment (adjust to your setup) ─────────────────────────────
# module load anaconda3
# conda activate snakemake_env
# source activate snakemake_env

# ─── working directory ────────────────────────────────────────────────────────
cd /home/usr/shenmj/2026-PhD_project
mkdir -p tmp/logs

echo "=============================="
echo "Job ID   : $SLURM_JOB_ID"
echo "Node     : $SLURMD_NODENAME"
echo "Start    : $(date)"
echo "=============================="

snakemake \
    -s /home/usr/shenmj/2026-PhD_project/script/snakemake.PhD2026.bamdam.all_batches.smk \
    -j 40 \
    --executor slurm \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going

echo "=============================="
echo "End : $(date)"
echo "=============================="
