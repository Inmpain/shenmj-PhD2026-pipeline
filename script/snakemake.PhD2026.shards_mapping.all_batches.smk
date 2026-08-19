# ─────────────────────────────────────────────────────────────────────────────
# YWLab PhD 2026 Project – Shards Mapping Pipeline (All Batches)
# Based on: snakemake.loess_batch1.shards_mapping.ver1.smk       (multi-lane)
#           snakemake_shards_mapping.ver1.250815.smk              (Lajia)
# Adapted for: shenmj | Date: 2026-06
#
# Input (from QC pipeline): {sample}.bbduk.lowcomp_filtered.fq
#
# Pipeline:
#   bowtie2_GTDB_mapping → bowtie2_shard_mapping (per shard, per type)
#
# Shard counts:
#   fungi:                8   (part_001 – part_008)
#   invertebrate:        100  (part_001 – part_100)
#   vertebrate-other:     80  (part_001 – part_080)
#   vertebrate-mammalian: 43  (part_001 – part_043)
#   plant:                80  (part_001 – part_080)
#   others:                2  (part_001 – part_002)
#
# Reference paths:
#   Multi-lane batches : /home/data/ref20250728/
#   Lajia              : /home/database/ref20250728/euk_ncbi_25Jul/
# ─────────────────────────────────────────────────────────────────────────────

import os

# ─────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────
RAW_BASE       = "/home/usr/shenmj/2026-PhD_project/00_raw_data"
PROCESSED_BASE = "/home/usr/shenmj/2026-PhD_project/01_processed_data"
BAM_BASE       = "/home/usr/shenmj/2026-PhD_project/03_bam"

# Reference databases
GTDB_REF_ML          = "/home/database/ref20250728/GTDB/GTDB.family_rep.bowtie2"
GTDB_REF_LAJIA       = "/home/database/ref20250728/GTDB/GTDB.family_rep.bowtie2"
SHARD_REF_BASE_ML    = "/home/database/ref20250728/euk_ncbi_25Jul"
SHARD_REF_BASE_LAJIA = "/home/database/ref20250728/euk_ncbi_25Jul"

# ─────────────────────────────────────────────────────────
# BATCH / SHARD CONFIGURATION
# ─────────────────────────────────────────────────────────
BATCH_LANES = {
    "GansuQinghai_samples_from202508": [f"L{i:03d}" for i in range(1, 9)],
    "GansuQinghai_samples_from202408": [f"L{i:03d}" for i in range(1, 5)],
    "Ningxia_samples"                : [f"L{i:03d}" for i in range(1, 5)],
}

SHARD_COUNTS = {
    "fungi"               : 8,
    "invertebrate"        : 100,
    "vertebrate-other"    : 80,
    "vertebrate-mammalian": 43,
    "plant"               : 80,
    "others"              : 2,
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
ML_MAPPING_TARGETS = []
for (batch, sn, snum) in VALID_ML_SAMPLES:
    for type_, n in SHARD_COUNTS.items():
        for num in [f"{i:03d}" for i in range(1, n + 1)]:
            ML_MAPPING_TARGETS.append(
                f"{BAM_BASE}/{batch}/{sn}_{snum}/{sn}_{snum}.{type_}.part_{num}.bam.finished"
            )

LAJIA_MAPPING_TARGETS = []
for sample in LAJIA_SAMPLES:
    for type_, n in SHARD_COUNTS.items():
        for num in [f"{i:03d}" for i in range(1, n + 1)]:
            LAJIA_MAPPING_TARGETS.append(
                f"{BAM_BASE}/Lajia_sites/{sample}/{sample}.{type_}.part_{num}.bam.finished"
            )

# ─────────────────────────────────────────────────────────
# ONSTART: environment / software check
# ─────────────────────────────────────────────────────────
onstart:
    import subprocess, sys
    print("\n" + "="*60)
    print("[ENV CHECK] Verifying required tools (mapping)...")
    _required = ["bowtie2", "samtools", "seqkit"]
    _missing  = []
    for _tool in _required:
        _r = subprocess.run(["which", _tool], capture_output=True)
        if _r.returncode == 0:
            _v = subprocess.run([_tool, "--version"], capture_output=True, text=True)
            _ver_line = (_v.stdout + _v.stderr).strip().split("\n")[0][:80]
            print(f"  ✓  {_tool:12s}  {_ver_line}")
        else:
            print(f"  ✗  {_tool:12s}  NOT FOUND")
            _missing.append(_tool)
    if _missing:
        sys.exit(f"\n[ENV CHECK] FATAL – tools not found: {_missing}\n"
                 "Activate your conda/module environment and retry.\n")
    os.makedirs(BAM_BASE, exist_ok=True)
    print("[ENV CHECK] All checks passed.\n" + "="*60 + "\n")

# ─────────────────────────────────────────────────────────
# WILDCARD CONSTRAINTS
# ─────────────────────────────────────────────────────────
wildcard_constraints:
    batch     = r"[A-Za-z0-9_]+",
    sn        = r"[^_/]+",
    snum      = r"S\d+",
    sample    = r"[^/]+",
    type      = r"fungi|invertebrate|vertebrate-other|vertebrate-mammalian|plant|others",
    shard_num = r"\d{3}"


# ═══════════════════════════════════════════════════════════════════════
# DEFAULT TARGET
# ═══════════════════════════════════════════════════════════════════════
rule all:
    input:
        ML_MAPPING_TARGETS + LAJIA_MAPPING_TARGETS


# ═══════════════════════════════════════════════════════════════════════
# MULTI-LANE MAPPING RULES
# Applies to: GansuQinghai_202508, GansuQinghai_202408, Ningxia
# ═══════════════════════════════════════════════════════════════════════

rule bowtie2_GTDB_mapping_ml:
    input:
        fq = PROCESSED_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.bbduk.lowcomp_filtered.fq"
    output:
        bam              = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.GTDB_family_rep.bam",
        mapped_flag      = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.GTDB_family_rep.bam.finished",
        reads_list       = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.GTDB_family_rep.mapped.reads",
        GTDB_filtered_fq = PROCESSED_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.non_GTDB.fq"
    params:
        out_dir            = BAM_BASE + "/{batch}/{sn}_{snum}/",
        ref_index_basename = GTDB_REF_ML
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        mkdir -p {params.out_dir}
        bowtie2 --threads {threads} \
        -x {params.ref_index_basename} -U {input.fq} -S {output.bam} \
        -k 100 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal \
        2> {output.bam}.log
        samtools view {output.bam}|cut -f1|uniq > {output.reads_list}
        seqkit grep -v -f {output.reads_list} {input.fq} > {output.GTDB_filtered_fq}
        touch {output.mapped_flag}
        """

rule bowtie2_shard_mapping_ml:
    input:
        fq = PROCESSED_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.non_GTDB.fq"
    output:
        bam         = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.{type}.part_{shard_num}.bam",
        mapped_flag = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.{type}.part_{shard_num}.bam.finished"
    params:
        ref_index_basename = SHARD_REF_BASE_ML + "/{type}.filtered.fa.split/{type}.filtered.part_{shard_num}.bowtie2"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        bowtie2 --threads {threads} \
        -x {params.ref_index_basename} -U {input.fq} -S {output.bam} \
        -k 100 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal \
        2> {output.bam}.log
        touch {output.mapped_flag}
        """


# ═══════════════════════════════════════════════════════════════════════
# LAJIA MAPPING RULES
# No lane structure; different ref path base
# ═══════════════════════════════════════════════════════════════════════

rule bowtie2_GTDB_mapping_lajia:
    input:
        fq = PROCESSED_BASE + "/Lajia_sites/{sample}/{sample}.bbduk.lowcomp_filtered.fq"
    output:
        bam              = BAM_BASE + "/Lajia_sites/{sample}/{sample}.GTDB_family_rep.bam",
        mapped_flag      = BAM_BASE + "/Lajia_sites/{sample}/{sample}.GTDB_family_rep.bam.finished",
        reads_list       = BAM_BASE + "/Lajia_sites/{sample}/{sample}.GTDB_family_rep.mapped.reads",
        GTDB_filtered_fq = PROCESSED_BASE + "/Lajia_sites/{sample}/{sample}.non_GTDB.fq"
    params:
        out_dir            = BAM_BASE + "/Lajia_sites/{sample}/",
        ref_index_basename = GTDB_REF_LAJIA
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        mkdir -p {params.out_dir}
        bowtie2 --threads {threads} \
        -x {params.ref_index_basename} -U {input.fq} -S {output.bam} \
        -k 100 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal \
        2> {output.bam}.log
        samtools view {output.bam}|cut -f1|uniq > {output.reads_list}
        seqkit grep -v -f {output.reads_list} {input.fq} > {output.GTDB_filtered_fq}
        touch {output.mapped_flag}
        """

rule bowtie2_shard_mapping_lajia:
    input:
        fq = PROCESSED_BASE + "/Lajia_sites/{sample}/{sample}.non_GTDB.fq"
    output:
        bam         = BAM_BASE + "/Lajia_sites/{sample}/{sample}.{type}.part_{shard_num}.bam",
        mapped_flag = BAM_BASE + "/Lajia_sites/{sample}/{sample}.{type}.part_{shard_num}.bam.finished"
    params:
        ref_index_basename = SHARD_REF_BASE_LAJIA + "/{type}.filtered.fa.split/{type}.filtered.part_{shard_num}.bowtie2"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        bowtie2 --threads {threads} \
        -x {params.ref_index_basename} -U {input.fq} -S {output.bam} \
        -k 100 -L 22 -i S,1,1.15 --mp 1,1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.1 --no-unal \
        2> {output.bam}.log
        touch {output.mapped_flag}
        """
