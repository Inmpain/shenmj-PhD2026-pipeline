#!/bin/bash
# =============================================================================
# YWLab PhD 2026 – submit NewDB metaDMG as one SLURM job per sample
#
# 遍历 03_bam 下所有样品目录，为每个还没有
# {sample}.newdb_cph_euk.metaDMG.aggregate.done 的样品投一个 metaDMG_newdb.slurm.sh
# 作业。可重复执行（已完成的自动跳过；已排队的按 jobname 去重）。
#
# 用法: bash submit_newdb_metaDMG_loop.sh
# 说明: 绕开 Snakemake 的 stat 轮询问题，逐样品独立作业。
# =============================================================================
set -euo pipefail

BAM_BASE=/home/usr/shenmj/2026-PhD_project/03_bam
LOG_DIR=/home/usr/shenmj/2026-PhD_project/tmp/logs
SCRIPT=/home/usr/shenmj/2026-PhD_project/script/metaDMG_newdb.slurm.sh

mkdir -p "$LOG_DIR"

submitted=0
skipped=0
no_sorted=0

for dir in $(ls -d "$BAM_BASE"/*/* 2>/dev/null); do
  batch=$(basename "$(dirname "$dir")")
  sample=$(basename "$dir")

  # 该样品必须已完成 sort（有 sorted.bam 且有 flag）
  if [ ! -f "$dir/$sample.newdb_cph_euk.merged.sorted.bam.flag" ]; then
    no_sorted=$((no_sorted+1))
    continue
  fi

  METADMG_BASE=/home/usr/shenmj/2026-PhD_project/04_metaDMG
  agg_flag="$METADMG_BASE/$batch/$sample/$sample.newdb_cph_euk.metaDMG.aggregate.done"

  if [ -f "$agg_flag" ]; then
    skipped=$((skipped+1))
    continue
  fi

  # 防重复提交：已排队/运行中的也跳过
  jobname="mdb_$(echo "$sample" | tr -cd 'A-Za-z0-9_' | cut -c1-20)"
  if squeue -u "$USER" --name="$jobname" -h 2>/dev/null | grep -q .; then
    echo "SKIP (already queued): $sample"
    skipped=$((skipped+1))
    continue
  fi

  sbatch \
    --job-name="$jobname" \
    --output="$LOG_DIR/metaDMG_${sample}.%j.log" \
    --error="$LOG_DIR/metaDMG_${sample}.%j.err" \
    "$SCRIPT" "$batch" "$sample"
  echo "SUBMIT $batch/$sample"
  submitted=$((submitted+1))
done

echo "============================================"
echo "submitted=$submitted  already_done_or_queued=$skipped  no_sorted=$no_sorted"
echo "aggregate done now: $(find /home/usr/shenmj/2026-PhD_project/04_metaDMG -name '*.newdb_cph_euk.metaDMG.aggregate.done' 2>/dev/null | wc -l)"
echo "total              : $(ls -d "$BAM_BASE"/*/* 2>/dev/null | wc -l)"
echo "============================================"