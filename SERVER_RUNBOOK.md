# YWLab PhD 2026 (shenmj) – alpha server paths & runbook

> Machine: `login01` (alpha), user `shenmj`
> Project root: `/home/usr/shenmj/2026-PhD_project`
> Storage note: `03_bam` is a symlink to `itp` storage node, **not mounted on node05/node06**
> → all sbatch jobs must carry `#SBATCH --exclude=node05,node06`

## 1. Project layout

| Path | Content |
|---|---|
| `00_raw_data/` | raw fastq, organized by batch (GansuQinghai_202508/202408, Ningxia_samples, Lajia_sites) |
| `01_processed_data/{batch}/{sample}/` | QC outputs, `{sample}.bbduk.lowcomp_filtered.fq`, `{sample}.non_GTDB.fq` |
| `02_qc/` | fastqc / seqkit QC artifacts + `.done` flags |
| `03_bam/` | **itp symlink** – GTDB bam, cph_euk newdb per-shard bams, merged/sorted bams |
| `04_metaDMG/` | metaDMG lca/dfit/aggregate outputs |
| `05_bamdam/` | bamdam outputs |
| `tmp/` | snakemake logs, `metaDMG_temp/`, `taxonomy_CPH/all.merged.acc2taxid` |
| `script/` | actual smk + sbatch files (this repo) |

## 2. Reference databases (alpha `/home/database/ref20250728/`)

### GTDB (microbial filter)
- `/home/database/ref20250728/GTDB/GTDB.family_rep.bowtie2`

### Old euk database (EOL legacy, replaced by cph_euk)
- `/home/database/ref20250728/euk_ncbi_25Jul`
- per-type shards: `{type}.filtered.fa.split/{type}.filtered.part_{num}.bowtie2`
- types: fungi(8) invertebrate(100) vertebrate-other(80) vertebrate-mammalian(43) plant(80) others(2)

### New euk database `cph_euk` (CURRENT)
- `/home/database/ref20250728/cph_euk`
- shards: `wgs_eukaryota.{1..129}.fas.gz` (bowtie2 index prefix)
- mapped per sample as `{sample}.newdb_cph_euk.wgs_eukaryota.{shard}.bam`

### Taxonomy (CURRENT, CPH build)
- root: `/home/database/ref20250728/taxonomy_CPH`
- `ncbi/20250530/nodes.dmp`, `ncbi/20250530/names.dmp` (used by metaDMG lca/dfit)
- acc2tax files (all merged into `tmp/taxonomy_CPH/all.merged.acc2taxid` by onstart):
  - `wgs_eukaryota.acc2taxid`
  - `cph_euk.plastid.mito.corent.acc2taxid`
  - `core_nt.acc2taxid`
  - `refseq_mitochondrion.genomic.acc2taxid`
  - `refseq_plastid.genomic.acc2taxid`
- `hires-organelles-viruses-smags/20240313/` (upstream hires taxonomy for microbial prefilter)

### Old taxonomy (EOL)
- `/home/database/ref20250728/taxonomy/` → `nodes_250309.dmp`, `names_250309.dmp`, `euk.acc2taxid`

### Misc
- adapter list: `/home/usr/xuez/adapter_list.fa`

## 3. Software / env

- snakemake: conda env `snakemake_env` (`module load anaconda3` / `conda activate snakemake_env`)
- tools: bowtie2, samtools, seqkit, fastp, bbduk, fastqc, metaDMG-cpp, bamdam
- SLURM executor: `--executor slurm -j 40 --latency-wait 60 --rerun-incomplete --keep-going`

## 4. Pipeline stage chain (current newdb path)

```
QC (bbduk.lowcomp_filtered.fq)
  → bowtie2 GTDB filter            -> {sample}.non_GTDB.fq
  → bowtie2 cph_euk 129 shards     -> {sample}.newdb_cph_euk.wgs_eukaryota.{1..129}.bam
  → merge + name-sort (qlen<200)   -> {sample}.newdb_cph_euk.merged.sorted.bam
  → metaDMG (lca/dfit/aggregate)   -> {sample}.newdb_cph_euk.metaDMG.aggregate.done
  → bamdam (upstream, not yet newdb-adapted)
```

## 5. Scripts actually used (in `script/`)

| File | Purpose |
|---|---|
| `snakemake.PhD2026.newdb_cph_euk_mapping.smk` | GTDB + cph_euk 129-shard mapping |
| `snakemake.PhD2026.newdb_cph_euk_postmapping.all_batches.smk` | merge + name-sort of 129 shards |
| `snakemake.PhD2026.newdb_cph_euk_metaDMG.all_batches.smk` | metaDMG on newdb sorted bam |
| `submit_newdb_postmapping.sbatch` | sort run (exclude node05/06) |
| `submit_newdb_metaDMG.sbatch` | metaDMG run (exclude node05/06) |
| `snakemake.PhD2026.shards_mapping.all_batches.smk` | EOL legacy euk mapping (313 shards) |
| `snakemake.PhD2026.shards_postmapping_qc_merge.all_batches.smk` | EOL legacy merge+sort |
| `snakemake.PhD2026.metaDMG.all_batches.smk` | EOL legacy metaDMG |
| `snakemake.PhD2026.bamdam.all_batches.smk` | bamdam (upstream; STRAND_TYPE=ds) |
| `snakemake.PhD2026.readsqc.REPAIR.smk` | QC repair for S76 (corrupt lane) + NTC S94 (empty) |
| `submit_*.sbatch` | corresponding drivers |

## 6. Verify sort completed (no error/no bug)

```bash
cd /home/usr/shenmj/2026-PhD_project

# 1) driver job exited clean
sacct -u shenmj --starttime today | grep -E "newdb_sort|COMPLETED|FAILED"
grep -iE "error|failed|exception" tmp/logs/snakemake_newdb_sort_*.err tmp/logs/snakemake_newdb_sort_*.log || echo "no errors in driver logs"

# 2) flag count == sample count
find 03_bam -name "*.newdb_cph_euk.merged.sorted.bam.flag" | wc -l

# 3) every sample dir has all 4 outputs (missing = not complete)
find 03_bam -name "*.newdb_cph_euk.merged.sorted.bam" | wc -l
find 03_bam -name "*.newdb_cph_euk.merged.sorted.bam.unique_mapped_reads" | wc -l

# 4) integrity: every sorted bam passes quickcheck
for f in $(find 03_bam -name "*.newdb_cph_euk.merged.sorted.bam"); do
  samtools quickcheck -v "$f" || echo "FAIL $f"
done
echo "quickcheck done"

# 5) per-sample job success inside snakemake log
grep -c "Finished job" tmp/logs/snakemake_newdb_sort_*.log
grep -E "Error in rule|Some jobs failed" tmp/logs/snakemake_newdb_sort_*.log || echo "no failed rules"

# 6) spot check header/reads of one bam
samtools view -H 03_bam/GansuQinghai_samples_from202508/*/*.newdb_cph_euk.merged.sorted.bam | head -3
```

## 7. Note on itp + node05/06

- `03_bam` lives on `itp`; node05/node06 do not mount `itp`
- even though `03_bam` is a symlink back to the original path on login nodes,
  compute jobs landing on node05/node06 would fail reading BAMs
- mitigation: `--exclude=node05,node06` on every sbatch driver
  (Snakemake SLURM executor submits per-job; the exclude is inherited by all child jobs)
