# ─────────────────────────────────────────────────────────────────────────────
# YWLab PhD 2026 Project – metaDMG Pipeline (All Batches)
# Based on: snakemake.loess_batch1.metaDMG.250817.smk  (multi-lane)
#           snakemake_metaDMG.250817.smk                (Lajia)
# Adapted for: shenmj | Date: 2026-06
#
# Input : {sample}.merged.sorted.bam  (from postmapping pipeline)
# Output: {sample}.metaDMG.aggregate.*
#
# Rule chain per sample:
#   metaDMG_dmg  ──┐
#   metaDMG_lca  ──┤
#   metaDMG_dfit ──┤
#                  └→ metaDMG_aggregate
#
# Reference paths:
#   Multi-lane batches : /home/data/ref20250728/
#   Lajia              : /home/database/ref20250728/taxonomy/
# ─────────────────────────────────────────────────────────────────────────────

import os

# ─────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────
RAW_BASE       = "/home/usr/shenmj/2026-PhD_project/00_raw_data"
BAM_BASE       = "/home/usr/shenmj/2026-PhD_project/03_bam"
METADMG_BASE   = "/home/usr/shenmj/2026-PhD_project/04_metaDMG"
METADMG_TEMP   = "/home/usr/shenmj/2026-PhD_project/tmp/metaDMG_temp"

# Reference databases
NODES_ML    = "/home/database/ref20250728/taxonomy/nodes_250309.dmp"
NAMES_ML    = "/home/database/ref20250728/taxonomy/names_250309.dmp"
ACC2TAX_ML  = "/home/database/ref20250728/taxonomy/euk.acc2taxid"

NODES_LAJIA   = "/home/database/ref20250728/taxonomy/nodes_250309.dmp"
NAMES_LAJIA   = "/home/database/ref20250728/taxonomy/names_250309.dmp"
ACC2TAX_LAJIA = "/home/database/ref20250728/taxonomy/euk.acc2taxid"

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
ML_METADMG_TARGETS = [
    f"{METADMG_BASE}/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.aggregate.done"
    for (batch, sn, snum) in VALID_ML_SAMPLES
]
LAJIA_METADMG_TARGETS = [
    f"{METADMG_BASE}/Lajia_sites/{sample}/{sample}.metaDMG.aggregate.done"
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
        ML_METADMG_TARGETS + LAJIA_METADMG_TARGETS


# ═══════════════════════════════════════════════════════════════════════
# MULTI-LANE metaDMG RULES
# ═══════════════════════════════════════════════════════════════════════

rule metaDMG_dmg_ml:
    input:
        merge_sort_bam = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam"
    output:
        flag = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.dmg.done"
    params:
        out_dir       = METADMG_BASE + "/{batch}/{sn}_{snum}/",
        output_prefix = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.dmg"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        mkdir -p {params.out_dir}
        metaDMG-cpp getdamage --run_mode 0 --print_length 15 --min_length 30 --out_prefix {params.output_prefix} {input.merge_sort_bam}
        touch {output.flag}
        """

rule metaDMG_lca_ml:
    input:
        merge_sort_bam = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam",
        nodes          = NODES_ML,
        names          = NAMES_ML,
        acc2tax        = ACC2TAX_ML,
    output:
        bdamage  = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.lca.bdamage.gz",
        lca_stat = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.lca.stat.gz",
        flag     = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.lca.done"
    params:
        temp_folder   = METADMG_TEMP + "/{batch}/{sn}_{snum}/",
        output_prefix = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.lca"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        # reallyDump 1 does not work
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
        names   = NAMES_ML,
    output:
        dfit = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.dfit.dfit.gz",
        flag = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.dfit.done"
    params:
        output_prefix = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.dfit"
    resources:
        mem_mb = 40000,
        nodes  = 1,
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
        names    = NAMES_ML,
    output:
        flag = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.aggregate.done"
    params:
        output_prefix = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.metaDMG.aggregate"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        # --dfit causes bug now!
        """
        metaDMG-cpp aggregate {input.bdamage} --dfit {input.dfit} --lcastat {input.lca_stat} --nodes {input.nodes} --names {input.names} --out {params.output_prefix}
        touch {output.flag}
        """

# Optional: ngsLCA (not connected to default target – run manually if needed)
rule ngsLCA_lca_ml:
    input:
        merge_sort_bam = BAM_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.merged.sorted.bam",
        nodes          = NODES_ML,
        names          = NAMES_ML,
        acc2tax        = ACC2TAX_ML,
    output:
        flag = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.ngsLCA.lca.done"
    params:
        output_prefix = METADMG_BASE + "/{batch}/{sn}_{snum}/{sn}_{snum}.ngsLCA.lca"
    shell:
        """
        ngsLCA -simscorelow 0.95 -simscorehigh 1.0 -fix_ncbi 0 -names {input.names} \
        -nodes {input.nodes} -acc2tax {input.acc2tax} -bam {input.merge_sort_bam} -outnames {params.output_prefix}
        """


# ═══════════════════════════════════════════════════════════════════════
# LAJIA metaDMG RULES
# ═══════════════════════════════════════════════════════════════════════

rule metaDMG_dmg_lajia:
    input:
        merge_sort_bam = BAM_BASE + "/Lajia_sites/{sample}/{sample}.merged.sorted.bam"
    output:
        flag = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.dmg.done"
    params:
        out_dir       = METADMG_BASE + "/Lajia_sites/{sample}/",
        output_prefix = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.dmg"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 20
    shell:
        """
        mkdir -p {params.out_dir}
        metaDMG-cpp getdamage --run_mode 0 --print_length 15 --min_length 30 --out_prefix {params.output_prefix} {input.merge_sort_bam}
        touch {output.flag}
        """

rule metaDMG_lca_lajia:
    input:
        merge_sort_bam = BAM_BASE + "/Lajia_sites/{sample}/{sample}.merged.sorted.bam",
        nodes          = NODES_LAJIA,
        names          = NAMES_LAJIA,
        acc2tax        = ACC2TAX_LAJIA,
    output:
        bdamage  = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.lca.bdamage.gz",
        lca_stat = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.lca.stat.gz",
        flag     = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.lca.done"
    params:
        temp_folder   = METADMG_TEMP + "/Lajia_sites/{sample}/",
        output_prefix = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.lca"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        # reallyDump 1 does not work
        """
        mkdir -p {params.temp_folder}
        metaDMG-cpp lca --threads {threads} --bam {input.merge_sort_bam} --nodes {input.nodes} --names {input.names} \
        --acc2tax {input.acc2tax} --fix_ncbi 0 --how_many 15 --sim_score_low 0.95 --weight_type 0 --lca_rank genus --temp {params.temp_folder} \
         --out_prefix {params.output_prefix}
        touch {output.flag}
        """

rule metaDMG_dfit_lajia:
    input:
        bdamage = rules.metaDMG_lca_lajia.output.bdamage,
        nodes   = NODES_LAJIA,
        names   = NAMES_LAJIA,
    output:
        dfit = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.dfit.dfit.gz",
        flag = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.dfit.done"
    params:
        output_prefix = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.dfit"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        """
        metaDMG-cpp dfit {input.bdamage} --threads {threads} --nodes {input.nodes} --names {input.names} \
        --nopt 5 --showfits 2 --seed 42 --out_prefix {params.output_prefix}
        touch {output.flag}
        """

rule metaDMG_aggregate_lajia:
    input:
        rules.metaDMG_dmg_lajia.output.flag,
        rules.metaDMG_dfit_lajia.output.flag,
        lca_stat = rules.metaDMG_lca_lajia.output.lca_stat,
        bdamage  = rules.metaDMG_lca_lajia.output.bdamage,
        dfit     = rules.metaDMG_dfit_lajia.output.dfit,
        nodes    = NODES_LAJIA,
        names    = NAMES_LAJIA,
    output:
        flag = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.aggregate.done"
    params:
        output_prefix = METADMG_BASE + "/Lajia_sites/{sample}/{sample}.metaDMG.aggregate"
    resources:
        mem_mb = 40000,
        nodes  = 1,
    threads: 40
    shell:
        # --dfit causes bug now!
        """
        metaDMG-cpp aggregate {input.bdamage} --dfit {input.dfit} --lcastat {input.lca_stat} --nodes {input.nodes} --names {input.names} --out {params.output_prefix}
        touch {output.flag}
        """
