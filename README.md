# shenmj-PhD2026-pipeline

Snakemake pipelines for YWLab PhD 2026 project (shenmj), running on the alpha cluster (login01).

## Status (2026-08)

- **cph_euk mapping DONE** (GTDB filter + 129 wgs_eukaryota shards) — output in `03_bam/{batch}/{sample}/`
- **postmapping (merge + name-sort) READY** — `submit_newdb_postmapping.sbatch`
- **metaDMG (newdb) READY** — `submit_newdb_metaDMG.sbatch`
- bamdam / legacy EOL (`euk_ncbi_25Jul`) / QC repair scripts also included for reference

## Layout

| Path | What |
|---|---|
| `script/*.smk` | snakemake workflows |
| `script/submit_*.sbatch` | SLURM driver scripts |
| `SERVER_RUNBOOK.md` | alpha server paths, database locations, verification commands |

## New-database pipeline (current)

```
QC -> GTDB filter -> cph_euk 129 shards -> merge+sort(qlen<200) -> metaDMG -> bamdam
```

Scripts actually used:
- `snakemake.PhD2026.newdb_cph_euk_mapping.smk`
- `snakemake.PhD2026.newdb_cph_euk_postmapping.all_batches.smk`
- `snakemake.PhD2026.newdb_cph_euk_metaDMG.all_batches.smk`
- `submit_newdb_postmapping.sbatch` / `submit_newdb_metaDMG.sbatch`

## Deploy to server

Public repo → direct `curl` download, e.g.:

```bash
mkdir -p /home/usr/shenmj/2026-PhD_project/script
cd /home/usr/shenmj/2026-PhD_project/script
curl -O https://raw.githubusercontent.com/Inmpain/shenmj-PhD2026-pipeline/main/script/snakemake.PhD2026.newdb_cph_euk_postmapping.all_batches.smk
curl -O https://raw.githubusercontent.com/Inmpain/shenmj-PhD2026-pipeline/main/script/snakemake.PhD2026.newdb_cph_euk_metaDMG.all_batches.smk
curl -O https://raw.githubusercontent.com/Inmpain/shenmj-PhD2026-pipeline/main/script/submit_newdb_postmapping.sbatch
curl -O https://raw.githubusercontent.com/Inmpain/shenmj-PhD2026-pipeline/main/script/submit_newdb_metaDMG.sbatch
chmod +x submit_newdb_*.sbatch
```

See `SERVER_RUNBOOK.md` for full path map + run verification.

## Notes

- `03_bam` is a symlink to the `itp` storage node, which is NOT mounted on node05/node06.
  All sbatch drivers therefore use `--exclude=node05,node06`.
- Params mirror the original group pipeline (Zhe Xue, 2025); see header comments in each `.smk`.
