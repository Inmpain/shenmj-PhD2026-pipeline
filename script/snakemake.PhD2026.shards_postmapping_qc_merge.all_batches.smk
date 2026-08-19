# ─────────────────────────────────────────────────────────────────────────────
# YWLab PhD 2026 Project – Post-mapping Merge & Sort Pipeline (All Batches)
# Based on: snakemake.loess_batch1.shards_postmapping_qc_merge.smk (multi-lane)
#           snakemake_shards_postmapping_qc_merge.smk               (Lajia)
# Adapted for: shenmj | Date: 2026-06
#
# Input : all per-shard BAMs from shards_mapping pipeline
# Output: {sample}.merged.sorted.bam  (reads <200 bp, for metaDMG)
#         {sample}.merged.sorted.bam.unique_mapped_reads
#
# Run AFTER snakemake.PhD2026.shards_mapping.all_batches.smk is complete.
# ─────────────────────────────────────────────────────────────────────────────

import os

# ─────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────
RAW_BASE = "/home/usr/shenmj/2026-PhD_project/00_raw_data"
BAM_BASE = "/home/usr/shenmj/2026-PhD_project/03_bam"

# ─────────────────────────────────────────────────────────
# BATCH CONFIGURATION
# ─────────────────────────────────────────────────────────
BATCH_LANES = {
    "GansuQinghai_samples_from202508": [f"L{i:03d}" for i in range(1, 9)],
    "GansuQinghai_samples_from202408": [f"L{i:03d}" for i in range(1, 5)],
    "Ningxia_samples"                : [f"L{i:03d}" for i in range(1, 5)],
}

# ─────────────────────────────────────────────────────────
# SAMPLE PARSING – multi-lane batches
# ─────────────────────────────────────────────────────────
VALID_ML_SAMPLES = {}

for _batch, _expected_lanes in BATCH_LANES.items():
    _raw_dir = os.path.join(RAW_BASE, _batch)
    if not os.path.isdir(_raw_dir):
        print(f"[WARNING] {_raw_dir} not found – skipping batch {_batch}")
        continue
    _per_sample = {}
    for _fname in os.listdir(_raw_dir):
        if not _fname.endswith("_R1_001.fastq.gz"):
            continue
        _parts = _fname.split("_")
        if len(_parts) < 5:
            continue
        _sn, _snum, _lane = _parts[0], _parts[1], _parts[2]
        _per_sample.setdefault((_sn, _snum), set()).add(_lane)
    _exp_set = set(_expected_lanes)
    for (_sn, _snum), _found in _per_sample.items():
        if _found == _exp_set:
            VALID_ML_SAMPLES[(_batch, _sn, _snum)] = sorted(_found)
        else:
            _miss  = sorted(_exp_set - _found)
            _extra = sorted(_found  - _exp_set)
            print(f"[WARNING] {_batch}/{_sn}_{_snum}: "
                  f"missing lanes {_miss}, unexpected {_extra} – SKIPPING")

# ─────────────────────────────────────────────────────────
# EXCLUDE NTC AND ZERO-READ SAMPLES FROM ALL DOWNSTREAM ANALYSIS
# LibNTC_S94: all 8 collapsed files are 0 bytes (expected NTC result).
# Bowtie2 / metaDMG / bamdam cannot process empty reads – excluded permanently.
# ─────────────────────────────────────────────────────────
_NTC_EXCLUDE = {
    ("GansuQinghai_samples_from202508", "LV7008875565-LibNTC25090803-LibNTC", "S94"),
}
_before = len(VALID_ML_SAMPLES)
VALID_ML_SAMPLES = {k: v for k, v in VALID_ML_SAMPLES.items() if k not in _NTC_EXCLUDE}
if _before - len(VALID_ML_SAMPLES) > 0:
    print(f"[INFO] Excluded NTC/zero-read samples: {_before - len(VALID_ML_SAMPLES)} "
          f"({[k[1]+'_'+k[2] for k in _NTC_EXCLUDE]})")

# ─────────────────────────────────────────────────────────
# SAMPLE PARSING – Lajia (no lanes, .fq.gz)
# ─────────────────────────────────────────────────────────
LAJIA_SAMPLES = []
_lajia_dir = os.path.join(RAW_BASE, "Lajia_sites")
if os.path.isdir(_lajia_dir):
    for _fname in os.listdir(_lajia_dir):
        if _fname.endswith("_R1.fq.gz"):
            LAJIA_SAMPLES.append(_fname[:-len("_R1.fq.gz")])
else:
    print(f"[WARNING] {_lajia_dir} not found – Lajia batch will be skipped")

print(f"[INFO] Multi-lane samples: {len(VALID_ML_SAMPLES)}")
print(f"[INFO] Lajia samples:      {len(LAJIA_SAMPLES)}")

# ─────────────────────────────────────────────────────────
# PER-SAMPLE BAM LISTS (wildcards escaped with double braces)
# ─────────────────────────────────────────────────────────
# Multi-lane
fungi_mapped_bams_ml = expand(
    BAM_BASE + "/{{batch}}/{{sn}}_{{snum}}/{{sn}}_{{snum}}.fungi.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 9)])
invert_mapped_bams_ml = expand(
    BAM_BASE + "/{{batch}}/{{sn}}_{{snum}}/{{sn}}_{{snum}}.invertebrate.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 101)])
vert_other_mapped_bams_ml = expand(
    BAM_BASE + "/{{batch}}/{{sn}}_{{snum}}/{{sn}}_{{snum}}.vertebrate-other.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 81)])
vert_mammal_mapped_bams_ml = expand(
    BAM_BASE + "/{{batch}}/{{sn}}_{{snum}}/{{sn}}_{{snum}}.vertebrate-mammalian.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 44)])
plant_mapped_bams_ml = expand(
    BAM_BASE + "/{{batch}}/{{sn}}_{{snum}}/{{sn}}_{{snum}}.plant.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 81)])
other_mapped_bams_ml = expand(
    BAM_BASE + "/{{batch}}/{{sn}}_{{snum}}/{{sn}}_{{snum}}.others.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 3)])

# Lajia
fungi_mapped_bams_lajia = expand(
    BAM_BASE + "/Lajia_sites/{{sample}}/{{sample}}.fungi.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 9)])
invert_mapped_bams_lajia = expand(
    BAM_BASE + "/Lajia_sites/{{sample}}/{{sample}}.invertebrate.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 101)])
vert_other_mapped_bams_lajia = expand(
    BAM_BASE + "/Lajia_sites/{{sample}}/{{sample}}.vertebrate-other.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 81)])
vert_mammal_mapped_bams_lajia = expand(
    BAM_BASE + "/Lajia_sites/{{sample}}/{{sample}}.vertebrate-mammalian.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 44)])
plant_mapped_bams_lajia = expand(
    BAM_BASE + "/Lajia_sites/{{sample}}/{{sample}}.plant.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 81)])
other_mapped_bams_lajia = expand(
    BAM_BASE + "/Lajia_sites/{{sample}}/{{sample}}.others.part_{num}.bam",
    num=[f"{i:03d}" for i in range(1, 3)])

# ─────────────────────────────────────────────────────────
# TARGET FILES
# ─────────────────────────────────────────────────────────
ML_MERGE_TARGETS = [
    f"{BAM_BASE}/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam.flag"
    for (batch, sn, snum) in VALID_ML_SAMPLES
]
LAJIA_MERGE_TARGETS = [
    f"{BAM_BASE}/Lajia_sites/{sample}/{sample}.merged.sorted.bam.flag"
    for sample in LAJIA_SAMPLES
]

# ─────────────────────────────────────────────────────────
# WILDCARD CONSTRAINTS
# ─────────────────────────────────────────────────────────
wildcard_constraints:
    batch  = r"[A-Za-z0-9_]+",
    sn     = r"[^_/]+",
    snum   = r"S\d+",
    sample = r"[^/]+"


# ═══════════════════════════════════════════════════════════════════════
# DEFAULT TARGET
# ═══════════════════════════════════════════════════════════════════════
rule all:
    input:
        ML_MERGE_TARGETS + LAJIA_MERGE_TARGETS


# ═══════════════════════════════════════════════════════════════════════
# MULTI-LANE MERGE & SORT
# ═══════════════════════════════════════════════════════════════════════

rule merge_sort_ml:
    input:
        bams = (fungi_mapped_bams_ml + invert_mapped_bams_ml +
                vert_other_mapped_bams_ml + vert_mammal_mapped_bams_ml +
                plant_mapped_bams_ml + other_mapped_bams_ml)
    output:
        merged_bam                = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.bam",
        merge_sort_bam            = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam",
        merge_sort_bam_uniq_reads = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam.unique_mapped_reads",
        merge_sort_bam_flag       = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam.flag",
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        samtools merge --threads {threads} -n -c -p -o {output.merged_bam} {input.bams}
        samtools sort -@ {threads} -n {output.merged_bam} | samtools view -b -h -e 'qlen<200' > {output.merge_sort_bam}
        samtools view {output.merge_sort_bam} |cut -f1|uniq > {output.merge_sort_bam_uniq_reads}
        touch {output.merge_sort_bam_flag}
        """


# ═══════════════════════════════════════════════════════════════════════
# LAJIA MERGE & SORT
# ═══════════════════════════════════════════════════════════════════════

rule merge_sort_lajia:
    input:
        bams = (fungi_mapped_bams_lajia + invert_mapped_bams_lajia +
                vert_other_mapped_bams_lajia + vert_mammal_mapped_bams_lajia +
                plant_mapped_bams_lajia + other_mapped_bams_lajia)
    output:
        merged_bam                = BAM_BASE + "/Lajia_sites/{sample}/{sample}.merged.bam",
        merge_sort_bam            = BAM_BASE + "/Lajia_sites/{sample}/{sample}.merged.sorted.bam",
        merge_sort_bam_uniq_reads = BAM_BASE + "/Lajia_sites/{sample}/{sample}.merged.sorted.bam.unique_mapped_reads",
        merge_sort_bam_flag       = BAM_BASE + "/Lajia_sites/{sample}/{sample}.merged.sorted.bam.flag",
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        samtools merge --threads {threads} -n -c -p -o {output.merged_bam} {input.bams}
        samtools sort -@ {threads} -n {output.merged_bam} | samtools view -b -h -e 'qlen<200' > {output.merge_sort_bam}
        samtools view {output.merge_sort_bam} |cut -f1|uniq > {output.merge_sort_bam_uniq_reads}
        touch {output.merge_sort_bam_flag}
        """
