# ThePipeline_v3

Lightweight, modular pipeline for processing MTBC whole-genome sequencing data: read cleaning, taxonomic filtering, mapping, variant calling and adjudication, consensus generation, resistance calling and typing.

This repository contains:

- `ThePipeline3` — main CLI driver that dispatches to modules in `PipeModules/`.
- `PipeModules/` — modular Python components for each pipeline stage (FastClean, Kraken, Mapping, Calling, Coverage, Consensus, Resistance, Typing, Distances, Clusters, MultiQC, Organize, plus helpers Repository/History/Version).
- `data/` — configuration and reference data (paths, configs, WHO catalog, annotations and vendored `pairsnp` library).
- `Programs/` — bundled program binaries and datasets (note: large/redistributable files have been removed; see `heavy_files_to_be_downloaded.txt`).
- `scripts/` — helper scripts, including `download_heavy_files.sh` to (optionally) re-acquire large external artifacts.
- `ThePipeline3_manual.md` — comprehensive manual auto-generated from source (parameters, file flows, examples).

Quick links

- Manual: `ThePipeline3_manual.md`
- Heavy files list: `heavy_files_to_be_downloaded.txt`
- Download helper: `scripts/download_heavy_files.sh`

Status

This repository is a working copy of ThePipeline3. Some very large external files (GATK jars, Minos Singularity image, and a snpEff data binary) were removed from the workspace to keep the repository small and reproducible. Use the download helper script or follow the instructions in `heavy_files_to_be_downloaded.txt` to obtain them before running the full pipeline.

Requirements

- Python 3.7
- Python packages: pandas, PyVCF (vcf), and other standard libraries. See `requirements.txt` for an environment snapshot.
- External tools (expected available under `Programs/` or on PATH): bwa, samtools, picard, fastp, seqtk, pigz, kraken, bedtools (genomeCoverageBed), gatk (Mutect2), VarScan (jar), Minos (Singularity image), snpEff, snp-sites, qualimap, MultiQC. Many of these are bundled under `Programs/` but large files are omitted — see below.

Getting ready (minimal)

1) Install Python packages (recommended via conda or venv). Example (using pip):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2) Re-acquire large external artifacts (GATK, Minos, snpEff data):

Use the helper script for guidance. For manual download follow instructions in `heavy_files_to_be_downloaded.txt`.

Automated (when you have direct URLs and permission):

```bash
MINOS_URL="https://example.com/minos_v0.12.5.img" \
GATK_LOCAL_URL="https://example.com/gatk-package-4.2.5.0-local.jar" \
GATK_SPARK_URL="https://example.com/gatk-package-4.2.5.0-spark.jar" \
SNIPEFF_BIN_URL="https://example.com/snpEffectPredictor.bin" \
./scripts/download_heavy_files.sh
```

Note: GATK distributions require license acceptance; you may need to download them manually.

Usage examples (dry-run commands are in the manual)

- Clean reads with fastp:

```bash
python3 ThePipeline3 fastclean -f sample_R1.fastq.gz sample_R2.fastq.gz -p sample -t 4
```

- Taxonomic filtering with Kraken:

```bash
python3 ThePipeline3 kraken -f sample.P1.clean.fastq.gz sample.P2.clean.fastq.gz --paired -p sample --classify --filter --report
```

- Mapping and variant calling:

```bash
python3 ThePipeline3 mapping -f sample.P1.filtered.fastq.gz sample.P2.filtered.fastq.gz -p sample -t 8 -c
python3 ThePipeline3 calling -p sample -t 4
```

Where to find more information

- Full parameter reference, inputs/outputs and examples: see `ThePipeline3_manual.md`.
- Program paths and data defaults are read from `data/Paths/programs_path` and `data/Paths/data_path` — edit these files if you install software in non-default locations.

Development notes

- Python scripts target Python 3.7 (shebangs use `#!/usr/bin/env python3.7`). Adjust the interpreter if necessary.
- The repo includes a vendored `pairsnp` implementation in `data/libs/pairsnp-python/` to avoid external package resolution issues.

Contributing

Open an issue or submit a PR for bugfixes and improvements. When adding new program versions, update `data/Configs/software_versions.txt` so history files reflect the new versions.

License

See `LICENSE` in the repository root.

Acknowledgements

This README and `ThePipeline3_manual.md` were generated programmatically by inspecting `ThePipeline3` and the `PipeModules/` folder for parameter names, defaults and I/O behaviour. Keep both documents in sync with code changes.
