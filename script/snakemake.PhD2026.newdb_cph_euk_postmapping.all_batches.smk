# ─────────────────────────────────────────────────────────────────────────────
# YWLab PhD 2026 Project – NewDB cph_euk Post-mapping Merge & Sort (All Batches)
# Based on:
#   snakemake.PhD2026.shards_postmapping_qc_merge.all_batches.smk  (old 313 shards)
#   new_single_multi/step4.euk.mapping.smk                         (robust merge+sort)
# Adapted for: shenmj | Date: 2026-08
#
# Input : 129 per-shard BAMs from snakemake.PhD2026.newdb_cph_euk_mapping.smk
#         {sample}.newdb_cph_euk.wgs_eukaryota.{shard}.bam  (shard 1-129)
# Output: {sample}.newdb_cph_euk.merged.bam                (name-merged, intermediate)
#         {sample}.newdb_cph_euk.merged.sorted.bam         (name-sorted, qlen<200)
#         {sample}.newdb_cph_euk.merged.sorted.bam.unique_mapped_reads
#
# Storage: BAM_BASE is symlinked to itp (not mounted on node05/node06)
#          -> submit via --exclude=node05,node06 (see submit_newdb_postmapping.sbatch)
#
# Run AFTER snakemake.PhD2026.newdb_cph_euk_mapping.smk is complete.
# GTDB outputs are NOT merged here (already filtered to non_GTDB.fq before cph_euk).
# ─────────────────────────────────────────────────────────────────────────────

import os
import glob

shell.executable("/bin/bash")

# ─────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────
PROJECT_BASE   = "/home/usr/shenmj/2026-PhD_project"
PROCESSED_BASE = os.path.join(PROJECT_BASE, "01_processed_data")
BAM_BASE       = os.path.join(PROJECT_BASE, "03_bam")

NEWDB_LABEL    = "newdb_cph_euk"
NEWDB_REF_BASE = "/home/database/ref20250728/cph_euk"  # reference only, not used in sort

# ─────────────────────────────────────────────────────────
# DATABASE SHARDS (must match newdb_cph_euk_mapping.smk)
# ─────────────────────────────────────────────────────────
NEWDB_SHARDS = [str(i) for i in range(1, 130)]  # 1-129 => 129 BAMs per sample

# ─────────────────────────────────────────────────────────
# SAMPLE DISCOVERY
# Mirrors snakemake.PhD2026.newdb_cph_euk_mapping.smk logic:
#   glob PROCESSED_BASE/*/*/*.bbduk.lowcomp_filtered.fq
#   (not RAW_BASE lane parsing, not per-batch hardcode)
# ─────────────────────────────────────────────────────────
QC_SUFFIX = ".bbduk.lowcomp_filtered.fq"

EXCLUDE_SAMPLES = {
    ("GansuQinghai_samples_from202508", "LV7008875565-LibNTC25090803-LibNTC_S94"),
}

SAMPLES = []

for fq in sorted(glob.glob(os.path.join(PROCESSED_BASE, "*", "*", f"*{QC_SUFFIX}"))):
    sample_dir = os.path.dirname(fq)
    sample = os.path.basename(sample_dir)
    batch = os.path.basename(os.path.dirname(sample_dir))
    expected_name = f"{sample}{QC_SUFFIX}"
    if os.path.basename(fq) != expected_name:
        continue
    if (batch, sample) in EXCLUDE_SAMPLES:
        continue
    SAMPLES.append((batch, sample))

SAMPLES = sorted(set(SAMPLES))

print(f"[INFO] NewDB postmapping – discovered {len(SAMPLES)} samples (129 shards each)")
for _b, _s in SAMPLES[:5]:
    print(f"  {_b}/{_s}")
if len(SAMPLES) > 5:
    print(f"  ... and {len(SAMPLES)-5} more")

# ─────────────────────────────────────────────────────────
# ONSTART
# ─────────────────────────────────────────────────────────
onstart:
    import subprocess, sys
    print("\n" + "="*60)
    print("[ENV CHECK] Verifying required tools (newdb postmapping)...")
    _required = ["samtools"]
    _missing = []
    for _tool in _required:
        _r = subprocess.run(["which", _tool], capture_output=True)
        if _r.returncode == 0:
            _v = subprocess.run([_tool, "--version"], capture_output=True, text=True)
            _ver = (_v.stdout + _v.stderr).strip().split("\n")[0][:80]
            print(f"  ✓  {_tool:12s}  {_ver}")
        else:
            print(f"  ✗  {_tool:12s}  NOT FOUND")
            _missing.append(_tool)
    if _missing:
        sys.exit(f"\n[ENV CHECK] FATAL – tools not found: {_missing}\n")
    os.makedirs(BAM_BASE, exist_ok=True)
    print(f"[ENV CHECK] BAM_BASE (itp symlink): {BAM_BASE}")
    print(f"[ENV CHECK] Expected BAMs per sample: {len(NEWDB_SHARDS)}")
    print("="*60 + "\n")

# ─────────────────────────────────────────────────────────
# HELPERS – expand per-sample BAM lists (wildcard-safe)
# ─────────────────────────────────────────────────────────
def _newdb_bams(wildcards):
    # called with {batch}/{sample} wildcards
    return expand(
        os.path.join(BAM_BASE, "{{batch}}", "{{sample}}", "{{sample}}." + NEWDB_LABEL + ".wgs_eukaryota.{shard}.bam"),
        shard=NEWDB_SHARDS
    )

def _newdb_flags(wildcards):
    return expand(
        os.path.join(BAM_BASE, "{{batch}}", "{{sample}}", "{{sample}}." + NEWDB_LABEL + ".wgs_eukaryota.{shard}.bam.finished"),
        shard=NEWDB_SHARDS
    )

# Used for rule input declaration (Snakemake expand with double braces for deferred wildcards)
newdb_bams_template = expand(
    BAM_BASE + "/{{batch}}/{{sample}}/{{sample}}." + NEWDB_LABEL + ".wgs_eukaryota.{shard}.bam",
    shard=NEWDB_SHARDS
)
newdb_flags_template = expand(
    BAM_BASE + "/{{batch}}/{{sample}}/{{sample}}." + NEWDB_LABEL + ".wgs_eukaryota.{shard}.bam.finished",
    shard=NEWDB_SHARDS
)

# ─────────────────────────────────────────────────────────
# TARGET FILES
# ─────────────────────────────────────────────────────────
MERGE_TARGETS = [
    os.path.join(BAM_BASE, batch, sample, f"{sample}.{NEWDB_LABEL}.merged.sorted.bam.flag")
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
        MERGE_TARGETS


# ═══════════════════════════════════════════════════════════════════════
# MERGE & SORT (NEWDB cph_euk)
# Robust pattern from new_single_multi/step4.euk.mapping.smk:
#   - BAM list file + expected count check
#   - samtools merge -b list | sort -n  in a pipe
#   - samtools view -e 'qlen<200' for aDNA length filter (kept from old postmapping)
#   - quickcheck + atomic mv via temporary file + trap
# ═══════════════════════════════════════════════════════════════════════
rule merge_sort_newdb:
    input:
        bams     = newdb_bams_template,
        finished = newdb_flags_template
    output:
        merged_bam  = os.path.join(BAM_BASE, "{batch}", "{sample}", "{sample}." + NEWDB_LABEL + ".merged.bam"),
        sorted_bam  = os.path.join(BAM_BASE, "{batch}", "{sample}", "{sample}." + NEWDB_LABEL + ".merged.sorted.bam"),
        uniq_reads  = os.path.join(BAM_BASE, "{batch}", "{sample}", "{sample}." + NEWDB_LABEL + ".merged.sorted.bam.unique_mapped_reads"),
        flag        = os.path.join(BAM_BASE, "{batch}", "{sample}", "{sample}." + NEWDB_LABEL + ".merged.sorted.bam.flag")
    log:
        os.path.join(BAM_BASE, "{batch}", "{sample}", "{sample}." + NEWDB_LABEL + ".merged.sorted.bam.log")
    params:
        out_dir       = os.path.join(BAM_BASE, "{batch}", "{sample}"),
        bam_list      = os.path.join(BAM_BASE, "{batch}", "{sample}", f".{{sample}}.{NEWDB_LABEL}.bam_list.txt"),
        tmp_merged    = os.path.join(BAM_BASE, "{batch}", "{sample}", f".{{sample}}.{NEWDB_LABEL}.merged.bam.tmp"),
        tmp_sorted    = os.path.join(BAM_BASE, "{batch}", "{sample}", f".{{sample}}.{NEWDB_LABEL}.merged.sorted.bam.tmp"),
        sort_tmp_prefix = os.path.join(BAM_BASE, "{batch}", "{sample}", f".{{sample}}.{NEWDB_LABEL}.sort_tmp"),
        expected_bams = len(NEWDB_SHARDS),
        merge_threads = 2,
        sort_threads  = 11,
        view_threads  = 2,
        sort_mem      = "2G"
    resources:
        mem_mb = 102400,  # 100 GiB, match step4 merge_and_name_sort
        nodes  = 1
    threads: 16
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.out_dir}
        mkdir -p "$(dirname {params.sort_tmp_prefix})"

        bam_list="{params.bam_list}"
        tmp_merged="{params.tmp_merged}"
        tmp_sorted="{params.tmp_sorted}"

        # Expand Snakemake's bam list into the file
        # NOTE: no :q on {input.bams} – for a list, :q quotes the WHOLE list as one
        # string (bam_list becomes 1 line and the 129-count check fails). Use plain join.
        printf '%s\n' {input.bams} > "${{bam_list}}"

        actual=$(wc -l < "${{bam_list}}")
        if [[ "${{actual}}" -ne {params.expected_bams} ]]; then
            echo "[ERROR] {wildcards.sample}: expected {params.expected_bams} BAMs, found ${{actual}}" >&2
            cat "${{bam_list}}" >&2
            exit 1
        fi

        trap 'rm -f "${{bam_list}}" "${{tmp_merged}}" "${{tmp_sorted}}"' EXIT
        rm -f "${{tmp_merged}}" "${{tmp_sorted}}"

        echo "[INFO] {wildcards.batch}/{wildcards.sample}: merging {params.expected_bams} shards..." | tee {log:q}
        echo "[INFO] BAM_BASE (itp symlink): {BAM_BASE}" | tee -a {log:q}

        # Validate inputs
        samtools quickcheck -v $(cat "${{bam_list}}") 2>&1 | tee -a {log:q}

        # Merge (name-collated) -> sort -n -> length filter
        # samtools merge produces name-collated when inputs are name-sorted; we enforce -n in sort
        samtools merge \
            -@ {params.merge_threads} \
            -u \
            -c -p --no-PG \
            -b "${{bam_list}}" \
            -o - 2>> {log} \
        | samtools sort \
            -n \
            -@ {params.sort_threads} \
            -m {params.sort_mem} \
            -T {params.sort_tmp_prefix} \
            -o "${{tmp_merged}}" \
            - 2>> {log}

        samtools quickcheck -v "${{tmp_merged}}" 2>&1 | tee -a {log}
        mv "${{tmp_merged}}" {output.merged_bam}

        # Length filter qlen<200 (aDNA) – keep header, filter reads
        samtools view \
            -@ {params.view_threads} \
            -b -h -e 'qlen<200' \
            {output.merged_bam} > "${{tmp_sorted}}" 2>> {log}

        samtools quickcheck -v "${{tmp_sorted}}" 2>&1 | tee -a {log}
        mv "${{tmp_sorted}}" {output.sorted_bam}

        # Unique mapped reads (for QC)
        samtools view {output.sorted_bam} | cut -f1 | uniq > {output.uniq_reads} 2>> {log}
        uniq_count=$(wc -l < {output.uniq_reads})
        echo "[INFO] {wildcards.sample}: unique_mapped_reads=${{uniq_count}}" | tee -a {log}

        touch {output.flag}
        trap - EXIT
        rm -f "${{bam_list}}"
        """
