# ─────────────────────────────────────────────────────────────────────────────
# YWLab PhD 2026 Project – Reads QC  *** REPAIR VERSION ***
#
# Fixes two problematic samples from GansuQinghai_samples_from202508:
#
#   1. LV7008875565-LibNTC25090803-LibNTC_S94  (NTC)
#      → All 8 lanes produced 0-byte collapsed files (expected: no real reads).
#        FastQC crashes on empty files. Fix: skip fastqc/seqkit when file is
#        empty, still touch all flags so the pipeline can complete.
#
#   2. LV7008875531-LV6000619333-YWL1-A7598_S76  (S76)
#      → L005 R2 raw file is 0 bytes (corrupt transfer). Other 7 lanes are OK.
#        Fix: skip L005, run merge_lane with L001-L004 + L006-L008 only.
#        Per-lane fastp/qc outputs for the 7 good lanes already exist.
#
# Run AFTER main QC pipeline (only these 2 samples will be processed).
# Usage:
#   Dry-run:  snakemake -s snakemake.PhD2026.readsqc.REPAIR.smk -j 40 --executor slurm --resources nodes=3 -n
#   Full run: snakemake -s snakemake.PhD2026.readsqc.REPAIR.smk -j 40 --executor slurm --resources nodes=3
# ─────────────────────────────────────────────────────────────────────────────

import os

# ─────────────────────────────────────────────────────────
# PATHS  (must match main pipeline)
# ─────────────────────────────────────────────────────────
ADAPTER_LIST   = "/home/usr/xuez/adapter_list.fa"
RAW_BASE       = "/home/usr/shenmj/2026-PhD_project/00_raw_data"
PROCESSED_BASE = "/home/usr/shenmj/2026-PhD_project/01_processed_data"
QC_BASE        = "/home/usr/shenmj/2026-PhD_project/02_qc"
STATS_DIR      = "/home/usr/shenmj/2026-PhD_project/tmp"

BATCH = "GansuQinghai_samples_from202508"

# ─────────────────────────────────────────────────────────
# REPAIR SAMPLE DEFINITIONS
# ─────────────────────────────────────────────────────────

# S76: 7 lanes only (L005 R2 is corrupt – skip it)
S76_SN   = "LV7008875531-LV6000619333-YWL1-A7598"
S76_SNUM = "S76"
S76_LANES = ["L001", "L002", "L003", "L004", "L006", "L007", "L008"]

# NTC: all 8 lanes, but with empty-file protection throughout
NTC_SN   = "LV7008875565-LibNTC25090803-LibNTC"
NTC_SNUM = "S94"
NTC_LANES = [f"L{i:03d}" for i in range(1, 9)]

print(f"[REPAIR] S76  lanes: {S76_LANES}")
print(f"[REPAIR] NTC  lanes: {NTC_LANES}")

# ─────────────────────────────────────────────────────────
# INPUT FUNCTIONS
# ─────────────────────────────────────────────────────────
def _s76_collapsed_fqs(wildcards):
    return expand(
        PROCESSED_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}_{{lane}}.fastp.collapsed.fq.gz",
        lane=S76_LANES)

def _s76_fastp_qc_flags(wildcards):
    return expand(
        QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}_{{lane}}.fastp.collapsed.fq.gz.qc.done",
        lane=S76_LANES)

def _s76_raw_qc_flags(wildcards):
    return expand(
        QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}_{{lane}}.raw.fq.fastqc.done",
        lane=S76_LANES)

def _ntc_collapsed_fqs(wildcards):
    return expand(
        PROCESSED_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}_{{lane}}.fastp.collapsed.fq.gz",
        lane=NTC_LANES)

def _ntc_fastp_qc_flags(wildcards):
    return expand(
        QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}_{{lane}}.fastp.collapsed.fq.gz.qc.done",
        lane=NTC_LANES)

def _ntc_raw_qc_flags(wildcards):
    return expand(
        QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}_{{lane}}.raw.fq.fastqc.done",
        lane=NTC_LANES)

# ─────────────────────────────────────────────────────────
# TARGETS
# ─────────────────────────────────────────────────────────
S76_TARGET = f"{QC_BASE}/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.allqc.done"
NTC_TARGET = f"{QC_BASE}/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.allqc.done"

# ═══════════════════════════════════════════════════════════════════════
# DEFAULT TARGET
# ═══════════════════════════════════════════════════════════════════════
rule all:
    input:
        [S76_TARGET, NTC_TARGET]


# ═══════════════════════════════════════════════════════════════════════
# S76 REPAIR RULES  (7-lane merge, skipping L005)
# Per-lane fastp + raw QC outputs for L001-L004,L006-L008 already exist.
# Only merge_lane onward needs to run.
# ═══════════════════════════════════════════════════════════════════════

rule merge_lane_s76:
    input:
        fqs       = _s76_collapsed_fqs,
        qc_flags  = _s76_fastp_qc_flags,
        raw_flags = _s76_raw_qc_flags,
    output:
        merged = PROCESSED_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.merged_lane.fq"
    resources:
        mem    = lambda w, attempt: f"{5 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        zcat {input.fqs} > {output.merged}
        """

rule qc_merge_lane_s76:
    input:
        merged = rules.merge_lane_s76.output.merged
    output:
        stats = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.merged_lane.fq.seqkit_stats",
        flag  = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.merged_lane.fq.qc.done"
    params:
        out_dir = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/"
    resources:
        mem    = lambda w, attempt: f"{20 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        fastqc --memory 10000 -f fastq -o {params.out_dir} {input.merged}
        seqkit stats -a {input.merged} > {output.stats}
        touch {output.flag}
        """

rule seqkit_derep_s76:
    input:
        merged   = rules.merge_lane_s76.output.merged,
        merge_qc = rules.qc_merge_lane_s76.output.flag,
    output:
        dereped = PROCESSED_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.seqkit.dereped.fq",
        flag    = PROCESSED_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.seqkit.dereped.fq.done"
    resources:
        mem    = lambda w, threads: f"{5 * threads} GiB",
        nodes  = 1,
    threads: 5
    shell:
        """
        seqkit rmdup --ignore-case --by-seq -o {output.dereped} {input.merged}
        touch {output.flag}
        """

rule qc_seqkit_derep_s76:
    input:
        dereped = rules.seqkit_derep_s76.output.dereped
    output:
        stats = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.seqkit.dereped.fq.seqkit_stats",
        flag  = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.seqkit.dereped.fq.fastqc.done"
    params:
        out_dir = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/"
    resources:
        mem    = lambda w, attempt: f"{10 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        fastqc --memory 10000 -f fastq -o {params.out_dir} {input.dereped}
        seqkit stats -a {input.dereped} > {output.stats}
        touch {output.flag}
        """

rule low_complexity_bbduk_s76:
    input:
        dereped = rules.seqkit_derep_s76.output.dereped
    output:
        filtered = PROCESSED_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.bbduk.lowcomp_filtered.fq",
        flag     = PROCESSED_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.bbduk.lowcomp_filtered.fq.done"
    resources:
        mem    = lambda w, threads: f"{5 * threads} GiB",
        nodes  = 1,
    shell:
        """
        bbduk.sh in={input.dereped} out={output.filtered} \
            maxns=25 minlen=30 entropy=0.7 entropywindow=30 entropyk=4
        touch {output.flag}
        """

rule qc_bbduk_s76:
    input:
        filtered = rules.low_complexity_bbduk_s76.output.filtered
    output:
        stats = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.bbduk.lowcomp_filtered.fq.seqkit_stats",
        flag  = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/{S76_SN}_{S76_SNUM}.bbduk.lowcomp_filtered.fq.fastqc.done"
    params:
        out_dir = QC_BASE + f"/{BATCH}/{S76_SN}_{S76_SNUM}/"
    resources:
        mem    = lambda w, attempt: f"{10 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        fastqc --memory 10000 -f fastq -o {params.out_dir} {input.filtered}
        seqkit stats -a {input.filtered} > {output.stats}
        touch {output.flag}
        """

rule qc_all_s76:
    input:
        rules.low_complexity_bbduk_s76.output.flag,
        rules.qc_bbduk_s76.output.flag,
        rules.qc_seqkit_derep_s76.output.flag,
    output:
        flag = S76_TARGET
    resources:
        mem   = lambda w, attempt: f"{1 * attempt} GiB",
        nodes = 1,
    shell:
        """
        touch {output.flag}
        """


# ═══════════════════════════════════════════════════════════════════════
# NTC REPAIR RULES  (empty-file protection throughout)
# fastp collapsed files already exist (0 bytes). Start from qc_fastp_merge.
# ═══════════════════════════════════════════════════════════════════════

rule qc_fastp_merge_ntc:
    # Run per-lane: skip fastqc/seqkit if collapsed file is empty
    input:
        collapsed = PROCESSED_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}_{{lane}}.fastp.collapsed.fq.gz"
    output:
        stats = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}_{{lane}}.fastp.collapsed.fq.gz.seqkit_stats",
        flag  = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}_{{lane}}.fastp.collapsed.fq.gz.qc.done"
    params:
        out_dir = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/"
    resources:
        mem    = lambda w, attempt: f"{5 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        if [ -s {input.collapsed} ]; then
            fastqc --memory 10000 -f fastq -o {params.out_dir} {input.collapsed}
            seqkit stats -a {input.collapsed} > {output.stats}
        else
            echo "WARN: empty collapsed file (NTC), skipping fastqc" > {output.stats}
        fi
        touch {output.flag}
        """

rule merge_lane_ntc:
    input:
        fqs       = _ntc_collapsed_fqs,
        qc_flags  = _ntc_fastp_qc_flags,
        raw_flags = _ntc_raw_qc_flags,
    output:
        merged = PROCESSED_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.merged_lane.fq"
    resources:
        mem    = lambda w, attempt: f"{5 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        zcat {input.fqs} > {output.merged} || true
        touch {output.merged}
        """

rule qc_merge_lane_ntc:
    input:
        merged = rules.merge_lane_ntc.output.merged
    output:
        stats = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.merged_lane.fq.seqkit_stats",
        flag  = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.merged_lane.fq.qc.done"
    params:
        out_dir = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/"
    resources:
        mem    = lambda w, attempt: f"{5 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        if [ -s {input.merged} ]; then
            fastqc --memory 10000 -f fastq -o {params.out_dir} {input.merged}
            seqkit stats -a {input.merged} > {output.stats}
        else
            echo "WARN: empty merged file (NTC), skipping fastqc" > {output.stats}
        fi
        touch {output.flag}
        """

rule seqkit_derep_ntc:
    input:
        merged   = rules.merge_lane_ntc.output.merged,
        merge_qc = rules.qc_merge_lane_ntc.output.flag,
    output:
        dereped = PROCESSED_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.seqkit.dereped.fq",
        flag    = PROCESSED_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.seqkit.dereped.fq.done"
    resources:
        mem    = lambda w, threads: f"{5 * threads} GiB",
        nodes  = 1,
    threads: 5
    shell:
        """
        if [ -s {input.merged} ]; then
            seqkit rmdup --ignore-case --by-seq -o {output.dereped} {input.merged}
        else
            touch {output.dereped}
        fi
        touch {output.flag}
        """

rule qc_seqkit_derep_ntc:
    input:
        dereped = rules.seqkit_derep_ntc.output.dereped
    output:
        stats = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.seqkit.dereped.fq.seqkit_stats",
        flag  = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.seqkit.dereped.fq.fastqc.done"
    params:
        out_dir = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/"
    resources:
        mem    = lambda w, attempt: f"{5 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        if [ -s {input.dereped} ]; then
            fastqc --memory 10000 -f fastq -o {params.out_dir} {input.dereped}
            seqkit stats -a {input.dereped} > {output.stats}
        else
            echo "WARN: empty dereped file (NTC), skipping fastqc" > {output.stats}
        fi
        touch {output.flag}
        """

rule low_complexity_bbduk_ntc:
    input:
        dereped = rules.seqkit_derep_ntc.output.dereped
    output:
        filtered = PROCESSED_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.bbduk.lowcomp_filtered.fq",
        flag     = PROCESSED_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.bbduk.lowcomp_filtered.fq.done"
    resources:
        mem    = lambda w, threads: f"{5 * threads} GiB",
        nodes  = 1,
    shell:
        """
        if [ -s {input.dereped} ]; then
            bbduk.sh in={input.dereped} out={output.filtered} \
                maxns=25 minlen=30 entropy=0.7 entropywindow=30 entropyk=4
        else
            touch {output.filtered}
        fi
        touch {output.flag}
        """

rule qc_bbduk_ntc:
    input:
        filtered = rules.low_complexity_bbduk_ntc.output.filtered
    output:
        stats = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.bbduk.lowcomp_filtered.fq.seqkit_stats",
        flag  = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/{NTC_SN}_{NTC_SNUM}.bbduk.lowcomp_filtered.fq.fastqc.done"
    params:
        out_dir = QC_BASE + f"/{BATCH}/{NTC_SN}_{NTC_SNUM}/"
    resources:
        mem    = lambda w, attempt: f"{5 * attempt} GiB",
        nodes  = 1,
    shell:
        """
        if [ -s {input.filtered} ]; then
            fastqc --memory 10000 -f fastq -o {params.out_dir} {input.filtered}
            seqkit stats -a {input.filtered} > {output.stats}
        else
            echo "WARN: empty bbduk output (NTC), skipping fastqc" > {output.stats}
        fi
        touch {output.flag}
        """

rule qc_all_ntc:
    input:
        rules.low_complexity_bbduk_ntc.output.flag,
        rules.qc_bbduk_ntc.output.flag,
        rules.qc_seqkit_derep_ntc.output.flag,
    output:
        flag = NTC_TARGET
    resources:
        mem   = lambda w, attempt: f"{1 * attempt} GiB",
        nodes = 1,
    shell:
        """
        touch {output.flag}
        """
