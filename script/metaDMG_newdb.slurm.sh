#!/bin/bash
# =============================================================================
# YWLab PhD 2026 – per-sample NewDB cph_euk metaDMG (single SLURM job)
# Usage: sbatch metaDMG_newdb.slurm.sh <batch> <sample>
#
# 四步链（与 snakemake.PhD2026.newdb_cph_euk_metaDMG.all_batches.smk 参数一致）：
#   getdamage -> lca -> dfit -> aggregate
# 幂等：若 {sample}.newdb_cph_euk.metaDMG.aggregate.done 已存在则直接退出。
# 绕开 Snakemake 在 itp/NFS 上的 stat 轮询（AttributeError: str has no
# is_storage / WorkflowError），逐样品独立作业更稳。
#
# 依赖：metaDMG-cpp, samtools；acc2tax 用已合并的 tmp/taxonomy_CPH/all.merged.acc2taxid
# 输出：04_metaDMG/{batch}/{sample}/ （本地路径，非软链）
# =============================================================================
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=40G
#SBATCH --exclude=node05,node06

set -euo pipefail

BATCH=$1
SAMPLE=$2

PROJECT_BASE=/home/usr/shenmj/2026-PhD_project
BAM_BASE=$PROJECT_BASE/03_bam
METADMG_BASE=$PROJECT_BASE/04_metaDMG
METADMG_TEMP=$PROJECT_BASE/tmp/metaDMG_temp

NODES=/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/nodes.dmp
NAMES=/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/names.dmp
ACC2TAX=$PROJECT_BASE/tmp/taxonomy_CPH/all.merged.acc2taxid

PREFIX=$SAMPLE.newdb_cph_euk
SORTED_BAM=$BAM_BASE/$BATCH/$SAMPLE/$PREFIX.merged.sorted.bam
OUT_DIR=$METADMG_BASE/$BATCH/$SAMPLE
TEMP_DIR=$METADMG_TEMP/$BATCH/$SAMPLE

AGG_FLAG=$OUT_DIR/$PREFIX.metaDMG.aggregate.done
DMG_FLAG=$OUT_DIR/$PREFIX.metaDMG.dmg.done
LCA_FLAG=$OUT_DIR/$PREFIX.metaDMG.lca.done
DFIT_FLAG=$OUT_DIR/$PREFIX.metaDMG.dfit.done

# 幂等：已完成则跳过
if [ -f "$AGG_FLAG" ]; then
  echo "[SKIP] $BATCH/$SAMPLE aggregate already done"
  exit 0
fi

# 输入校验
if [ ! -f "$SORTED_BAM" ]; then
  echo "[ERROR] $BATCH/$SAMPLE: sorted bam missing: $SORTED_BAM" >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$TEMP_DIR"

echo "[INFO] $BATCH/$SAMPLE: step 1/4 getdamage ..."

metaDMG-cpp getdamage \
    --run_mode 0 \
    --print_length 15 \
    --min_length 30 \
    --out_prefix "$OUT_DIR/$PREFIX.metaDMG.dmg" \
    "$SORTED_BAM"
touch "$DMG_FLAG"

echo "[INFO] $BATCH/$SAMPLE: step 2/4 lca ..."

metaDMG-cpp lca \
    --threads ${SLURM_CPUS_PER_TASK:-40} \
    --bam "$SORTED_BAM" \
    --nodes "$NODES" \
    --names "$NAMES" \
    --acc2tax "$ACC2TAX" \
    --fix_ncbi 0 \
    --how_many 15 \
    --sim_score_low 0.95 \
    --weight_type 0 \
    --lca_rank genus \
    --temp "$TEMP_DIR" \
    --out_prefix "$OUT_DIR/$PREFIX.metaDMG.lca"
touch "$LCA_FLAG"

echo "[INFO] $BATCH/$SAMPLE: step 3/4 dfit ..."

metaDMG-cpp dfit "$OUT_DIR/$PREFIX.metaDMG.lca.bdamage.gz" \
    --threads ${SLURM_CPUS_PER_TASK:-40} \
    --nodes "$NODES" \
    --names "$NAMES" \
    --nopt 5 \
    --showfits 2 \
    --seed 42 \
    --out_prefix "$OUT_DIR/$PREFIX.metaDMG.dfit"
touch "$DFIT_FLAG"

echo "[INFO] $BATCH/$SAMPLE: step 4/4 aggregate ..."

metaDMG-cpp aggregate "$OUT_DIR/$PREFIX.metaDMG.lca.bdamage.gz" \
    --dfit "$OUT_DIR/$PREFIX.metaDMG.dfit.dfit.gz" \
    --lcastat "$OUT_DIR/$PREFIX.metaDMG.lca.stat.gz" \
    --nodes "$NODES" \
    --names "$NAMES" \
    --out "$OUT_DIR/$PREFIX.metaDMG.aggregate"
touch "$AGG_FLAG"

echo "[DONE] $BATCH/$SAMPLE metaDMG complete"