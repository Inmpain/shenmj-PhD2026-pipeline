# ─────────────────────────────────────────────────────────────────────────────
# YWLab PhD 2026 Project – NewDB cph_euk metaDMG Pipeline (All Batches)
# Based on: snakemake.PhD2026.metaDMG.all_batches.smk  (old euk, euk.acc2taxid)
#           new_single_multi/step5.euk.metaDMG.smk      (wgs+refseq mito/plastid)
# Adapted for: shenmj | Date: 2026-08
#
# Input : {sample}.newdb_cph_euk.merged.sorted.bam  (from newdb postmapping)
# Output: {sample}.newdb_cph_euk.metaDMG.aggregate.done  (+ lca/dfit/dmg flags)
#
# Reference:
#   NODES/NAMES : /home/database/ref20250728/taxonomy_CPH/ncbi/20250530/nodes.dmp
#   ACC2TAX     : 全量合并 taxonomy_CPH 下 5 个 acc2taxid
#                 wgs_eukaryota + cph_euk.plastid.mito.corent
#                 + core_nt + refseq_mitochondrion.genomic + refseq_plastid.genomic
#                (header dedup via awk, see onstart check)
#                -> 合并产物: tmp/taxonomy_CPH/all.merged.acc2taxid
#
# Params: 完全沿用旧 metaDMG.all_batches.smk，不改阈值
#   dmg : getdamage --run_mode 0 --print_length 15 --min_length 30
#   lca : --fix_ncbi 0 --how_many 15 --sim_score_low 0.95 --weight_type 0 --lca_rank genus
#   dfit: --nopt 5 --showfits 2 --seed 42
# ─────────────────────────────────────────────────────────────────────────────

import os
import glob

# ─────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────
PROJECT_BASE   = "/home/usr/shenmj/2026-PhD_project"
PROCESSED_BASE = os.path.join(PROJECT_BASE, "01_processed_data")
BAM_BASE       = os.path.join(PROJECT_BASE, "03_bam")
METADMG_BASE   = os.path.join(PROJECT_BASE, "04_metaDMG")
METADMG_TEMP   = os.path.join(PROJECT_BASE, "tmp/metaDMG_temp")

# ---- Taxonomy (alpha / CPH) ----
# 你在 taxonomy_CPH 下 tree 所示：ncbi/20250530 为 wgs 专用
NODES_ML      = "/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/nodes.dmp"
NAMES_ML      = "/home/database/ref20250728/taxonomy_CPH/ncbi/20250530/names.dmp"

# ---- acc2tax ----
# 全部并入（wgs + cph_euk plastid/mito/corent + core_nt + refseq mito/plastid）
WGS_ACC2TAX              = "/home/database/ref20250728/taxonomy_CPH/wgs_eukaryota.acc2taxid"
CPH_PLASTID_MITO_ACC2TAX = "/home/database/ref20250728/taxonomy_CPH/cph_euk.plastid.mito.corent.acc2taxid"
CORE_NT_ACC2TAX          = "/home/database/ref20250728/taxonomy_CPH/core_nt.acc2taxid"
REFSEQ_MITO_ACC2TAX      = "/home/database/ref20250728/taxonomy_CPH/refseq_mitochondrion.genomic.acc2taxid"
REFSEQ_PLASTID_ACC2TAX   = "/home/database/ref20250728/taxonomy_CPH/refseq_plastid.genomic.acc2taxid"

ACC2TAX_SOURCES = [
    WGS_ACC2TAX,
    CPH_PLASTID_MITO_ACC2TAX,
    CORE_NT_ACC2TAX,
    REFSEQ_MITO_ACC2TAX,
    REFSEQ_PLASTID_ACC2TAX,
]

# 合并后供 lca 使用（自动去表头、去重）
# 放项目 tmp 下，保证 shenmj 有写权限（database 目录未必可写）
MERGED_ACC2TAX           = os.path.join(PROJECT_BASE, "tmp/taxonomy_CPH/all.merged.acc2taxid")

# lca rule 引用的 acc2tax
ACC2TAX_ML = MERGED_ACC2TAX

# 新库 BAM 前缀（必须与 postmapping 输出一致）
NEWDB_LABEL = "newdb_cph_euk"

# ─────────────────────────────────────────────────────────
# SAMPLE DISCOVERY – 与 newdb postmapping 完全一致 (glob, 非 lane)
# ─────────────────────────────────────────────────────────
QC_SUFFIX = ".bbduk.lowcomp_filtered.fq"
EXCLUDE_SAMPLES = {
    ("GansuQinghai_samples_from202508", "LV7008875565-LibNTC25090803-LibNTC_S94"),
}

SAMPLES = []
for fq in sorted(glob.glob(os.path.join(PROCESSED_BASE, "*", "*", f"*{QC_SUFFIX}"))):
    sd = os.path.dirname(fq)
    sm = os.path.basename(sd)
    bt = os.path.basename(os.path.dirname(sd))
    if os.path.basename(fq) != f"{sm}{QC_SUFFIX}":
        continue
    if (bt, sm) in EXCLUDE_SAMPLES:
        continue
    SAMPLES.append((bt, sm))
SAMPLES = sorted(set(SAMPLES))

print(f"[INFO] NewDB metaDMG – total {len(SAMPLES)} samples (Lajia batch handled via {{batch}}/{sample})")
print(f"[INFO] NODES: {NODES_ML}")
print(f"[INFO] MERGED_ACC2TAX: {MERGED_ACC2TAX}")
print(f"[INFO]   acc2tax sources ({len(ACC2TAX_SOURCES)}):")
for _p in ACC2TAX_SOURCES:
    print(f"[INFO]     - {_p}")

# ─────────────────────────────────────────────────────────
# ONSTART: 检查并按需合并 acc2tax (防多次表头)
# ─────────────────────────────────────────────────────────
onstart:
    import subprocess, sys, pathlib
    print("\n" + "="*60)
    print("[ENV CHECK] metaDMG (newdb_cph_euk) – taxonomy & tools")
    # tools
    for _tool in ["metaDMG-cpp", "samtools"]:
        _r = subprocess.run(["which", _tool], capture_output=True)
        if _r.returncode == 0:
            _v = subprocess.run([_tool, "--help"], capture_output=True, text=True)
            _ver = (_v.stdout + _v.stderr).strip().split("\n")[0][:80]
            print(f"  ✓  {_tool:12s}  {_ver}")
        else:
            print(f"  ✗  {_tool:12s}  NOT FOUND")
    # taxonomy files
    for _p in [NODES_ML, NAMES_ML] + ACC2TAX_SOURCES:
        if not pathlib.Path(_p).is_file():
            print(f"  ✗  missing: {_p}")
        else:
            print(f"  ✓  {_p}")
    # merged acc2tax – 若不存在或需更新，现场合并（去表头 + 去重提示）
    _merged = pathlib.Path(MERGED_ACC2TAX)
    _need_merge = not _merged.is_file()
    # 若任一源文件新于 merged，则提示重建
    if _merged.is_file():
        _mt = _merged.stat().st_mtime
        for _src in ACC2TAX_SOURCES:
            try:
                if pathlib.Path(_src).stat().st_mtime > _mt:
                    print(f"[WARN] source newer than merged: {_src} ->建议重建 MERGED_ACC2TAX")
            except: pass
    if _need_merge:
        print(f"[INFO] MERGED_ACC2TAX not found, building: {MERGED_ACC2TAX}")
        # 去表头：若首行含 'accession' 则跳过；否则全量合并
        # 用 awk 处理，避免 python 大文件内存问题
        os.makedirs(str(pathlib.Path(MERGED_ACC2TAX).parent), exist_ok=True)
        _src_list = " ".join(f'"{p}"' for p in ACC2TAX_SOURCES)
        _cmd = f"""awk 'FNR==1 && tolower($0) ~ /^accession/ {{next}} 1' {_src_list} > "{MERGED_ACC2TAX}.tmp" && mv "{MERGED_ACC2TAX}.tmp" "{MERGED_ACC2TAX}" && wc -l "{MERGED_ACC2TAX}" """
        print(f"[CMD] {_cmd}")
        _r = subprocess.run(_cmd, shell=True, capture_output=True, text=True)
        print(_r.stdout)
        if _r.stderr:
            print(_r.stderr)
        if _r.returncode != 0:
            sys.exit(f"[ENV CHECK] FATAL – failed to build MERGED_ACC2TAX")
        print(f"  ✓  built {MERGED_ACC2TAX}")
    else:
        print(f"  ✓  merged exists: {MERGED_ACC2TAX} ({_merged.stat().st_size/1024/1024:.1f} MB)")
    # 确保输出根
    os.makedirs(METADMG_BASE, exist_ok=True)
    os.makedirs(METADMG_TEMP, exist_ok=True)
    print("="*60 + "\n")

# ─────────────────────────────────────────────────────────
# TARGET FILES
# 统一用 {batch}/{sample} 通配符（Lajia 的 batch 即 "Lajia_sites"），
# 不再拆分 ML / Lajia 两套规则，避免 AmbiguousRuleException。
# ─────────────────────────────────────────────────────────
ALL_METADMG_TARGETS = [
    os.path.join(METADMG_BASE, batch, sample, f"{sample}.{NEWDB_LABEL}.metaDMG.aggregate.done")
    for batch, sample in SAMPLES
]

# ─────────────────────────────────────────────────────────
# WILDCARD CONSTRAINTS
# ─────────────────────────────────────────────────────────
wildcard_constraints:
    batch  = r"[^/]+",
    sample = r"[^/]+"

# ═══════════════════════════════════════════════════════════════════════
# DEFAULT TARGET
# ═══════════════════════════════════════════════════════════════════════
rule all:
    input:
        ALL_METADMG_TARGETS


# ═══════════════════════════════════════════════════════════════════════
# MULTI-LANE metaDMG RULES – 参数完全沿用 snakemake.PhD2026.metaDMG.all_batches.smk
# 输入 BAM 已改为 newdb_cph_euk.merged.sorted.bam
# 统一 {batch}/{sample} 通配符（Lajia 的 batch 即 "Lajia_sites"），不拆分 ML/Lajia 规则
# ═══════════════════════════════════════════════════════════════════════

rule metaDMG_dmg_ml:
    input:
        merge_sort_bam = BAM_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".merged.sorted.bam"
    output:
        flag = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.dmg.done"
    params:
        out_dir       = METADMG_BASE + "/{batch}/{sample}/",
        output_prefix = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.dmg"
    resources:
        mem_mb = 40000,
        nodes  = 1
    threads: 20
    shell:
        """
        mkdir -p {params.out_dir}
        metaDMG-cpp getdamage --run_mode 0 --print_length 15 --min_length 30 --out_prefix {params.output_prefix} {input.merge_sort_bam}
        touch {output.flag}
        """

rule metaDMG_lca_ml:
    input:
        merge_sort_bam = BAM_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".merged.sorted.bam",
        nodes          = NODES_ML,
        names          = NAMES_ML,
        acc2tax        = ACC2TAX_ML
    output:
        bdamage  = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.lca.bdamage.gz",
        lca_stat = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.lca.stat.gz",
        flag     = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.lca.done"
    params:
        temp_folder   = METADMG_TEMP + "/{batch}/{sample}/",
        output_prefix = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.lca"
    resources:
        mem_mb = 40000,
        nodes  = 1
    threads: 40
    shell:
        """
        mkdir -p {params.temp_folder}
        metaDMG-cpp lca --threads {threads} --bam {input.merge_sort_bam} --nodes {input.nodes} --names {input.names} \
        --acc2tax {input.acc2tax} --fix_ncbi 0 --how_many 15 --sim_score_low 0.95 --weight_type 0 --lca_rank genus --temp {params.temp_folder} \
         --out_prefix {params.output_prefix}
        touch {output.flag}
        """

rule metaDMG_dfit_ml:
    input:
        bdamage = rules.metaDMG_lca_ml.output.bdamage,
        nodes   = NODES_ML,
        names   = NAMES_ML
    output:
        dfit = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.dfit.dfit.gz",
        flag = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.dfit.done"
    params:
        output_prefix = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.dfit"
    resources:
        mem_mb = 40000,
        nodes  = 1
    threads: 40
    shell:
        """
        metaDMG-cpp dfit {input.bdamage} --threads {threads} --nodes {input.nodes} --names {input.names} \
        --nopt 5 --showfits 2 --seed 42 --out_prefix {params.output_prefix}
        touch {output.flag}
        """

rule metaDMG_aggregate_ml:
    input:
        rules.metaDMG_dmg_ml.output.flag,
        rules.metaDMG_dfit_ml.output.flag,
        lca_stat = rules.metaDMG_lca_ml.output.lca_stat,
        bdamage  = rules.metaDMG_lca_ml.output.bdamage,
        dfit     = rules.metaDMG_dfit_ml.output.dfit,
        nodes    = NODES_ML,
        names    = NAMES_ML
    output:
        flag = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.aggregate.done"
    params:
        output_prefix = METADMG_BASE + "/{batch}/{sample}/{sample}." + NEWDB_LABEL + ".metaDMG.aggregate"
    resources:
        mem_mb = 40000,
        nodes  = 1
    threads: 40
    shell:
        """
        metaDMG-cpp aggregate {input.bdamage} --dfit {input.dfit} --lcastat {input.lca_stat} --nodes {input.nodes} --names {input.names} --out {params.output_prefix}
        touch {output.flag}
        """
