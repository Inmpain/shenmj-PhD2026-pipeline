#!/bin/bash
#SBATCH --job-name=PhD26_newdb_metaDMG
#SBATCH --output=/home/usr/shenmj/2026-PhD_project/tmp/logs/snakemake_newdb_metaDMG_%j.log
#SBATCH --error=/home/usr/shenmj/2026-PhD_project/tmp/logs/snakemake_newdb_metaDMG_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# ─── BAM on itp symlink (03_bam); metaDMG reads via 03_bam path ─────────────
# ─── node05/node06 don't mount itp -> keep exclude for metaDMG reads too ────
#SBATCH --exclude=node05,node06

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
echo "Script   : snakemake.PhD2026.newdb_cph_euk_metaDMG.all_batches.smk"
echo "Acc2Tax  : merged (wgs_eukaryota + cph_euk.plastid.mito.corent)"
echo "=============================="

snakemake \
    -s /home/usr/shenmj/2026-PhD_project/script/snakemake.PhD2026.newdb_cph_euk_metaDMG.all_batches.smk \
    -j 40 \
    --executor slurm \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going

echo "=============================="
echo "End : $(date)"
echo "=============================="
