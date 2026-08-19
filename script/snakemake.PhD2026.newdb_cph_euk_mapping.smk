# =============================================================================
# PhD 2026 project — GTDB + cph_euk NewDB mapping
#
# Workflow:
#   sample.bbduk.lowcomp_filtered.fq
#       -> GTDB mapping
#       -> sample.non_GTDB.fq
#       -> cph_euk NewDB mapping (wgs_eukaryota.1–129)
#
# Existing complete GTDB outputs are reused automatically.
# =============================================================================

import os
import glob

shell.executable("/bin/bash")


# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------

PROJECT_BASE = "/home/usr/shenmj/2026-PhD_project"

PROCESSED_BASE = os.path.join(PROJECT_BASE, "01_processed_data")
BAM_BASE       = os.path.join(PROJECT_BASE, "03_bam")

GTDB_REF = "/home/database/ref20250728/GTDB/GTDB.family_rep.bowtie2"

NEWDB_REF_BASE = "/home/database/ref20250728/cph_euk"
NEWDB_LABEL    = "newdb_cph_euk"


# -----------------------------------------------------------------------------
# RESOURCE SETTINGS
# -----------------------------------------------------------------------------

GTDB_THREADS = 20
GTDB_MEM_MB  = 40000

NEWDB_BOWTIE2_THREADS  = 20
NEWDB_SAMTOOLS_THREADS = 4
NEWDB_THREADS          = NEWDB_BOWTIE2_THREADS + NEWDB_SAMTOOLS_THREADS

# 可根据首次大 shard 测试结果调整
NEWDB_MEM_MB = 200000


# -----------------------------------------------------------------------------
# DATABASE SHARDS
# -----------------------------------------------------------------------------

NEWDB_SHARDS = [str(i) for i in range(1, 130)]


# -----------------------------------------------------------------------------
# SAMPLE DISCOVERY
#
# 只从已经 QC 完成的样品级文件发现样品：
#   01_processed_data/{batch}/{sample}/{sample}.bbduk.lowcomp_filtered.fq
#
# 不扫描原始 fastq；不检查 lane；不区分 ML / Lajia。
# -----------------------------------------------------------------------------

QC_SUFFIX = ".bbduk.lowcomp_filtered.fq"

# 保留旧流程中明确排除的零 reads NTC
EXCLUDE_SAMPLES = {
    (
        "GansuQinghai_samples_from202508",
        "LV7008875565-LibNTC25090803-LibNTC_S94"
    ),
}

SAMPLES = []

for fq in sorted(
    glob.glob(os.path.join(PROCESSED_BASE, "*", "*", f"*{QC_SUFFIX}"))
):
    sample_dir = os.path.dirname(fq)
    sample = os.path.basename(sample_dir)
    batch = os.path.basename(os.path.dirname(sample_dir))

    expected_name = f"{sample}{QC_SUFFIX}"

    # 只接受标准目录和标准命名的样品级 QC 文件
    if os.path.basename(fq) != expected_name:
        continue

    if (batch, sample) in EXCLUDE_SAMPLES:
        continue

    SAMPLES.append((batch, sample))

SAMPLES = sorted(set(SAMPLES))


# -----------------------------------------------------------------------------
# FINAL TARGETS
# -----------------------------------------------------------------------------

NEWDB_TARGETS = [
    os.path.join(
        BAM_BASE,
        batch,
        sample,
        f"{sample}.{NEWDB_LABEL}.wgs_eukaryota.{shard}.bam.finished"
    )
    for batch, sample in SAMPLES
    for shard in NEWDB_SHARDS
]


# -----------------------------------------------------------------------------
# WILDCARDS
# -----------------------------------------------------------------------------

wildcard_constraints:
    batch       = r"[^/]+",
    sample      = r"[^/]+",
    newdb_shard = r"\d+"


# -----------------------------------------------------------------------------
# DEFAULT TARGET
# -----------------------------------------------------------------------------

rule all:
    input:
        NEWDB_TARGETS


# -----------------------------------------------------------------------------
# GTDB MAPPING
#
# 保持原脚本的 GTDB 参数、输出命名和 non_GTDB 生成逻辑不变。
# 若四个 output 均已存在，Snakemake 自动跳过此 rule。
# -----------------------------------------------------------------------------

rule bowtie2_GTDB_mapping:
    input:
        fq = os.path.join(
            PROCESSED_BASE,
            "{batch}",
            "{sample}",
            "{sample}.bbduk.lowcomp_filtered.fq"
        )

    output:
        bam = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}.GTDB_family_rep.bam"
        ),
        finished = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}.GTDB_family_rep.bam.finished"
        ),
        reads_list = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}.GTDB_family_rep.mapped.reads"
        ),
        non_gtdb_fq = os.path.join(
            PROCESSED_BASE,
            "{batch}",
            "{sample}",
            "{sample}.non_GTDB.fq"
        )

    params:
        out_dir = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}"
        ),
        ref_index_basename = GTDB_REF

    threads:
        GTDB_THREADS

    resources:
        mem_mb = GTDB_MEM_MB,
        nodes = 1

    log:
        os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}.GTDB_family_rep.bam.log"
        )

    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.out_dir}

        bowtie2 \
            --threads {threads} \
            -x {params.ref_index_basename} \
            -U {input.fq} \
            -S {output.bam} \
            -k 100 \
            -L 22 \
            -i S,1,1.15 \
            --mp 1,1 \
            --rdg 0,1 \
            --rfg 0,1 \
            --score-min L,0,-0.1 \
            --no-unal \
            2> {log}

        samtools view {output.bam} \
            | cut -f1 \
            | uniq \
            > {output.reads_list}

        seqkit grep \
            -v \
            -f {output.reads_list} \
            {input.fq} \
            > {output.non_gtdb_fq}

        touch {output.finished}
        """


# -----------------------------------------------------------------------------
# cph_euk NEW DATABASE MAPPING
#
# non_GTDB.fq 和 GTDB finished marker 都依赖上一条 GTDB rule。
# 因而：
#   已完成 GTDB 的样品 -> 直接进入 NewDB mapping
#   未完成 GTDB 的样品 -> 先跑 GTDB，再进入 NewDB mapping
# -----------------------------------------------------------------------------

rule bowtie2_newdb_mapping:
    input:
        fq = os.path.join(
            PROCESSED_BASE,
            "{batch}",
            "{sample}",
            "{sample}.non_GTDB.fq"
        ),
        gtdb_done = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}.GTDB_family_rep.bam.finished"
        )

    output:
        bam = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}." + NEWDB_LABEL
            + ".wgs_eukaryota.{newdb_shard}.bam"
        ),
        finished = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}." + NEWDB_LABEL
            + ".wgs_eukaryota.{newdb_shard}.bam.finished"
        )

    params:
        out_dir = os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}"
        ),
        ref_index_basename = os.path.join(
            NEWDB_REF_BASE,
            "wgs_eukaryota.{newdb_shard}.fas.gz"
        ),
        bowtie2_threads = NEWDB_BOWTIE2_THREADS,
        samtools_threads = NEWDB_SAMTOOLS_THREADS

    threads:
        NEWDB_THREADS

    resources:
        mem_mb = NEWDB_MEM_MB,
        nodes = 1

    log:
        os.path.join(
            BAM_BASE,
            "{batch}",
            "{sample}",
            "{sample}." + NEWDB_LABEL
            + ".wgs_eukaryota.{newdb_shard}.bam.log"
        )

    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.out_dir}

        bowtie2 \
            --threads {params.bowtie2_threads} \
            -x {params.ref_index_basename} \
            -U {input.fq} \
            -k 100 \
            -L 22 \
            -i S,1,1.15 \
            --mp 1,1 \
            --rdg 0,1 \
            --rfg 0,1 \
            --score-min L,0,-0.1 \
            --no-unal \
            2> {log} \
        | samtools view \
            -@ {params.samtools_threads} \
            -bh \
            -o {output.bam} \
            -

        samtools quickcheck -v {output.bam}

        touch {output.finished}
        """
