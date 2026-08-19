#!/bin/bash
#SBATCH --job-name=PhD26_newdb_sort
#SBATCH --output=/home/usr/shenmj/2026-PhD_project/tmp/logs/snakemake_newdb_sort_%j.log
#SBATCH --error=/home/usr/shenmj/2026-PhD_project/tmp/logs/snakemake_newdb_sort_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --exclude=node05,node06

# ─── BAM storage on itp (symlinked to 03_bam) is not mounted on node05/node06 ─
#     --exclude ensures merge jobs never land on those nodes

# ─── activate environment ──────────────────────────────────────────────────
# module load samtools
# conda activate snakemake_env

# ─── working directory ─────────────────────────────────────────────────────
cd /home/usr/shenmj/2026-PhD_project
mkdir -p tmp/logs

echo "=============================="
echo "Job ID   : $SLURM_JOB_ID"
echo "Node     : $SLURMD_NODENAME"
echo "Start    : $(date)"
echo "BAM_BASE : /home/usr/shenmj/2026-PhD_project/03_bam (itp symlink)"
echo "Script   : snakemake.PhD2026.newdb_cph_euk_postmapping.all_batches.smk"
echo "=============================="

snakemake \
    -s /home/usr/shenmj/2026-PhD_project/script/snakemake.PhD2026.newdb_cph_euk_postmapping.all_batches.smk \
    -j 40 \
    --executor slurm \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going

echo "=============================="
echo "End : $(date)"
echo "=============================="
