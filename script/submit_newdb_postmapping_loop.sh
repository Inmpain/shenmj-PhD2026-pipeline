#!/bin/bash
# =============================================================================
# YWLab PhD 2026 – submit NewDB postmapping (merge+sort) as one SLURM job per sample
#
# 遍历 03_bam 下所有样品目录，为每个还没有
# {sample}.newdb_cph_euk.merged.sorted.bam.flag 的样品投一个
# merge_sort_newdb.slurm.sh 作业。可重复执行（已完成的自动跳过）。
#
# 用法: bash submit_newdb_postmapping_loop.sh
# 说明: 绕开 Snakemake 在 itp/NFS 上的 stat 轮询（WorkflowError），逐样品独立作业更稳。
# =============================================================================
set -euo pipefail

BAM_BASE=/home/usr/shenmj/2026-PhD_project/03_bam
LOG_DIR=/home/usr/shenmj/2026-PhD_project/tmp/logs
SCRIPT=/home/usr/shenmj/2026-PhD_project/script/merge_sort_newdb.slurm.sh

mkdir -p "$LOG_DIR"

submitted=0
skipped=0
failed_dir=0

for dir in $(ls -d "$BAM_BASE"/*/* 2>/dev/null); do
  batch=$(basename "$(dirname "$dir")")
  sample=$(basename "$dir")

  # 目录本身要真的是样品（含 shard bam），排除日志/临时目录等
  if [ ! -f "$dir/$sample.newdb_cph_euk.wgs_eukaryota.1.bam" ]; then
    echo "SKIP (no shard bam): $dir"
    failed_dir=$((failed_dir+1))
    continue
  fi

  flag="$dir/$sample.newdb_cph_euk.merged.sorted.bam.flag"
  if [ -f "$flag" ]; then
    skipped=$((skipped+1))
    continue
  fi

  # 防止重复提交：已排队/运行中的也跳过
  if squeue -u "$USER" --name="ms_${sample:0:20}" -h 2>/dev/null | grep -q .; then
    echo "SKIP (already queued): $sample"
    skipped=$((skipped+1))
    continue
  fi

  jobname="ms_$(echo "$sample" | tr -cd 'A-Za-z0-9_' | cut -c1-24)"
  sbatch \
    --job-name="$jobname" \
    --output="$LOG_DIR/merge_${sample}.%j.log" \
    --error="$LOG_DIR/merge_${sample}.%j.err" \
    "$SCRIPT" "$batch" "$sample"
  echo "SUBMIT $batch/$sample"
  submitted=$((submitted+1))
done

echo "============================================"
echo "submitted=$submitted  already_done_or_queued=$skipped  no_shard_dir=$failed_dir"
echo "done now: $(ls "$BAM_BASE"/*/*/*.newdb_cph_euk.merged.sorted.bam.flag 2>/dev/null | wc -l)"
echo "total   : $(ls -d "$BAM_BASE"/*/*/*.newdb_cph_euk.wgs_eukaryota.1.bam 2>/dev/null | wc -l)"
echo "============================================"