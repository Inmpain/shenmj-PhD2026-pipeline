#!/bin/bash
# =============================================================================
# YWLab PhD 2026 – per-sample NewDB cph_euk merge + name-sort (single SLURM job)
# Usage: sbatch merge_sort_newdb.slurm.sh <batch> <sample>
#
# 与 snakemake.PhD2026.newdb_cph_euk_postmapping.all_batches.smk 的
# merge_sort_newdb 规则完全等价（同参数同输出命名），
# 但作为独立 slurm 作业运行，绕开 Snakemake 在 itp/NFS 上的持续 stat 轮询
# （WorkflowError: unable to obtain modification time ...）。
#
# 依赖：samtools
# 输出（均在 itp 软链 03_bam 下）：
#   {sample}.newdb_cph_euk.merged.bam
#   {sample}.newdb_cph_euk.merged.sorted.bam           (name-sorted, qlen<200)
#   {sample}.newdb_cph_euk.merged.sorted.bam.unique_mapped_reads
#   {sample}.newdb_cph_euk.merged.sorted.bam.flag
# =============================================================================
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --exclude=node05,node06

set -euo pipefail

BATCH=$1
SAMPLE=$2

BAM_BASE=/home/usr/shenmj/2026-PhD_project/03_bam
OUT_DIR=$BAM_BASE/$BATCH/$SAMPLE
PREFIX=$SAMPLE.newdb_cph_euk
EXPECTED=129

LOG=$OUT_DIR/$PREFIX.merged.sorted.bam.log
BAM_LIST=$OUT_DIR/.$PREFIX.bam_list.txt
TMP_MERGED=$OUT_DIR/.$PREFIX.merged.bam.tmp
TMP_SORTED=$OUT_DIR/.$PREFIX.merged.sorted.bam.tmp
SORT_TMP=$OUT_DIR/.$PREFIX.sort_tmp

MERGED_BAM=$OUT_DIR/$PREFIX.merged.bam
SORTED_BAM=$OUT_DIR/$PREFIX.merged.sorted.bam
UNIQ_READS=$OUT_DIR/$PREFIX.merged.sorted.bam.unique_mapped_reads
FLAG=$OUT_DIR/$PREFIX.merged.sorted.bam.flag

mkdir -p "$OUT_DIR" "$(dirname "$SORT_TMP")"

# ---- 收集实际存在的 shard bam（防缺 shard 时 brace-expansion 假计数） ----
> "$BAM_LIST"
for i in $(seq 1 "$EXPECTED"); do
  if [ -f "$OUT_DIR/$PREFIX.wgs_eukaryota.$i.bam" ]; then
    echo "$OUT_DIR/$PREFIX.wgs_eukaryota.$i.bam" >> "$BAM_LIST"
  fi
done

actual=$(wc -l < "$BAM_LIST")
if [ "$actual" -ne "$EXPECTED" ]; then
  echo "[ERROR] $SAMPLE: expected $EXPECTED BAMs, found $actual" >&2
  cat "$BAM_LIST" >&2
  exit 1
fi

trap 'rm -f "$BAM_LIST" "$TMP_MERGED" "$TMP_SORTED"' EXIT
rm -f "$TMP_MERGED" "$TMP_SORTED"

echo "[INFO] $BATCH/$SAMPLE: merging $EXPECTED shards..." | tee "$LOG"
echo "[INFO] BAM_BASE (itp symlink): $BAM_BASE" | tee -a "$LOG"

# ---- 输入完整性 ----
samtools quickcheck -v $(cat "$BAM_LIST") 2>&1 | tee -a "$LOG"

# ---- merge (name-collated) -> sort -n ----
samtools merge \
    -@ 2 \
    -u \
    -c -p --no-PG \
    -b "$BAM_LIST" \
    -o - 2>> "$LOG" \
| samtools sort \
    -n \
    -@ 11 \
    -m 2G \
    -T "$SORT_TMP" \
    -o "$TMP_MERGED" \
    - 2>> "$LOG"

samtools quickcheck -v "$TMP_MERGED" 2>&1 | tee -a "$LOG"
mv "$TMP_MERGED" "$MERGED_BAM"

# ---- aDNA 长度过滤 qlen<200 ----
samtools view \
    -@ 2 \
    -b -h -e 'qlen<200' \
    "$MERGED_BAM" > "$TMP_SORTED" 2>> "$LOG"

samtools quickcheck -v "$TMP_SORTED" 2>&1 | tee -a "$LOG"
mv "$TMP_SORTED" "$SORTED_BAM"

# ---- unique mapped reads ----
samtools view "$SORTED_BAM" | cut -f1 | uniq > "$UNIQ_READS" 2>> "$LOG"
uniq_count=$(wc -l < "$UNIQ_READS")
echo "[INFO] $SAMPLE: unique_mapped_reads=$uniq_count" | tee -a "$LOG"

touch "$FLAG"
trap - EXIT
rm -f "$BAM_LIST"