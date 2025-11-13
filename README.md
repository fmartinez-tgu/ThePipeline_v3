# ThePipeline_v3

Lightweight, modular pipeline for processing MTBC whole-genome sequencing data: read cleaning, taxonomic filtering, mapping, variant calling and adjudication, consensus generation, resistance calling and typing.

This repository contains:

- `ThePipeline3` — main CLI (Command Line Interface) driver that dispatches to modules in `PipeModules/`.
- `PipeModules/` — modular Python components for each pipeline stage (FastClean, Kraken, Mapping, Calling, Coverage, Consensus, Resistance, Typing, Distances, Clusters, MultiQC, Organize, plus helpers Repository/History/Version).
- `data/` — configuration and reference data (paths, configs, WHO catalog, annotations and vendored `pairsnp` library).
- `Programs/` — bundled program binaries and datasets (note: large/redistributable files have been removed; see `heavy_files_to_be_downloaded.txt`).
- `scripts/` — helper scripts, including `download_heavy_files.sh` to (optionally) re-acquire large external artifacts.
- `ThePipeline3_manual.md` — all you need to know about this Pipeline: parameters, file flows and examples.

**Requirements**

- Python 3.7
- Python packages: pandas, PyVCF (vcf), and other standard libraries. See `DEPENDENCIES.txt`.
- External tools (expected available under `Programs/` or on PATH): bwa, samtools, picard, fastp, seqtk, pigz, kraken, bedtools (genomeCoverageBed), gatk (Mutect2), VarScan, Minos (Singularity image), snpEff, snp-sites, qualimap, MultiQC. Many of these are bundled under `Programs/` but large files are omitted — see below.


**Usage examples (dry-run commands are in the manual)**

- Clean reads with fastp:

```bash
ThePipeline3 fastclean -f sample_R1.fastq.gz sample_R2.fastq.gz -p sample -t 4 -v
```

- Taxonomic filtering with Kraken:

```bash
ThePipeline3 kraken -f sample.P1.clean.fastq.gz sample.P2.clean.fastq.gz --paired --compressed -p sample --classify --filter --report -t 5
```

- Mapping:

```bash
ThePipeline3 mapping -f sample.P1.filtered.fastq.gz sample.P2.filtered.fastq.gz -p sample -t 8 -c
```

- Coverage filtering:
```bash
ThePipeline3 coverage -p sample -e <file format, e.g. cram> -f 
```

- Variant calling:
```bash
ThePipeline3 calling -p sample -t 4 -e <file format> 
```

**Where to find more information**

- Full parameter reference, inputs/outputs and examples: see `ThePipeline3_manual.md`.
- Program paths and data defaults are read from `data/Paths/programs_path` and `data/Paths/data_path` — edit these files if you install software in non-default locations.

**Development notes**

- Python scripts target Python 3.7. Adjust the interpreter if necessary.
- The repo includes a vendored `pairsnp` implementation in `data/libs/pairsnp-python/` to avoid external package resolution issues.

**Contributing**

Open an issue or submit a PR for bugfixes and improvements. When adding new program versions, update `data/Configs/software_versions.txt` so history files reflect the new versions.

**License**

See `LICENSE` in the repository root.
