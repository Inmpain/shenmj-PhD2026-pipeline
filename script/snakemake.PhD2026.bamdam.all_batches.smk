# ─────────────────────────────────────────────────────────────────────────────
# YWLab PhD 2026 Project – BamDam Pipeline (All Batches)
# Based on: snakemake.bamdam.smk  (lake_huron / single-lane)
# Adapted for: shenmj | Date: 2026-06
#
# Input : {sample}.merged.sorted.bam  (from postmapping pipeline)
#         {sample}.metaDMG.lca.stat.gz (from metaDMG pipeline)
# Output: {sample}.bamdam.compute.stats / .subs
#
# Rule chain per sample:
#   bamdam_shrink → bamdam_compute
#
# !!!! IMPORTANT: SET STRAND_TYPE BEFORE RUNNING !!!!
#   STRAND_TYPE = "ds"  for double-stranded libraries
#   STRAND_TYPE = "ss"  for single-stranded libraries
# ─────────────────────────────────────────────────────────────────────────────

import os

############################################
#######IMPORTANT!!!!!!!!!!!!!!!!!!!!########
############################################

##SET THIS TO CORRECT TYPE BEFORE RUNNING!!!
##OPTIONS ARE ss OR ds######################
STRAND_TYPE = "ds"

############################################
#######READ ABOVE###########################
############################################

# ─────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────
RAW_BASE     = "/home/usr/shenmj/2026-PhD_project/00_raw_data"
BAM_BASE     = "/home/usr/shenmj/2026-PhD_project/03_bam"
METADMG_BASE = "/home/usr/shenmj/2026-PhD_project/04_metaDMG"
BAMDAM_BASE  = "/home/usr/shenmj/2026-PhD_project/05_bamdam"

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
# TARGET FILES
# ─────────────────────────────────────────────────────────
ML_BAMDAM_TARGETS = [
    f"{BAMDAM_BASE}/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.compute.done"
    for (batch, sn, snum) in VALID_ML_SAMPLES
]
LAJIA_BAMDAM_TARGETS = [
    f"{BAMDAM_BASE}/Lajia_sites/{sample}/{sample}.bamdam.compute.done"
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
        ML_BAMDAM_TARGETS + LAJIA_BAMDAM_TARGETS


# ═══════════════════════════════════════════════════════════════════════
# MULTI-LANE BamDam RULES
# ═══════════════════════════════════════════════════════════════════════

rule bamdam_shrink_ml:
    input:
        merge_sort_bam = BAM_BASE     + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam",
        lca_stat       = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.lca.stat.gz",
    output:
        shrink_bam = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.shrink.bam",
        shrink_lca = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.shrink.lca",
        flag       = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.shrink.done"
    params:
        out_dir       = BAMDAM_BASE + "/{batch}/{sn}_{snum}/",
        output_prefix = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.shrink"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        """
        mkdir -p {params.out_dir}
        bamdam shrink --in_lca {input.lca_stat} --in_bam {input.merge_sort_bam} --out_lca {output.shrink_lca} --out_bam {output.shrink_bam} \
        --stranded {STRAND_TYPE} --mincount 3 --upto species --show_progress
        touch {output.flag}
        """

rule bamdam_compute_ml:
    input:
        shrink_bam = rules.bamdam_shrink_ml.output.shrink_bam,
        shrink_lca = rules.bamdam_shrink_ml.output.shrink_lca,
    output:
        stats = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.compute.stats",
        subs  = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.compute.subs",
        flag  = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.compute.done"
    params:
        output_prefix = BAMDAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bamdam.compute"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        """
        bamdam compute --in_lca {input.shrink_lca} --in_bam {input.shrink_bam} --out_tsv {output.stats} --out_subs {output.subs} \
        --stranded {STRAND_TYPE} --upto species --show_progress
        touch {output.flag}
        """


# ═══════════════════════════════════════════════════════════════════════
# LAJIA BamDam RULES
# ═══════════════════════════════════════════════════════════════════════

rule bamdam_shrink_lajia:
    input:
        merge_sort_bam = BAM_BASE     + "/Lajia_sites/{sample}/{sample}.merged.sorted.bam",
        lca_stat       = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.lca.stat.gz",
    output:
        shrink_bam = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.shrink.bam",
        shrink_lca = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.shrink.lca",
        flag       = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.shrink.done"
    params:
        out_dir       = BAMDAM_BASE + "/Lajia_sites/{sample}/",
        output_prefix = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.shrink"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        """
        mkdir -p {params.out_dir}
        bamdam shrink --in_lca {input.lca_stat} --in_bam {input.merge_sort_bam} --out_lca {output.shrink_lca} --out_bam {output.shrink_bam} \
        --stranded {STRAND_TYPE} --mincount 3 --upto species --show_progress
        touch {output.flag}
        """

rule bamdam_compute_lajia:
    input:
        shrink_bam = rules.bamdam_shrink_lajia.output.shrink_bam,
        shrink_lca = rules.bamdam_shrink_lajia.output.shrink_lca,
    output:
        stats = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.compute.stats",
        subs  = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.compute.subs",
        flag  = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.compute.done"
    params:
        output_prefix = BAMDAM_BASE + "/Lajia_sites/{sample}/{sample}.bamdam.compute"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        """
        bamdam compute --in_lca {input.shrink_lca} --in_bam {input.shrink_bam} --out_tsv {output.stats} --out_subs {output.subs} \
        --stranded {STRAND_TYPE} --upto species --show_progress
        touch {output.flag}
        """
