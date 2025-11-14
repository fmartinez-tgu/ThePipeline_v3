ThePipeline3 — Reference Manual
================================

This manual documents ThePipeline3 (Python-based), its submodules in `PipeModules/`, the command-line parameters exposed by `ThePipeline3` driver, the required inputs, intermediate files produced during execution and which are removed, and the final outputs.

Contents
--------
- [About](#about)
- [Short Dependencies](#short-dependencies)
- [Layout and configuration files imported by the pipeline](#layout-and-configuration-files-imported-by-the-pipeline)
- [How the driver works (ThePipeline3)](#how-the-driver-works-thepipeline3)
- [Per-module reference (parameters, inputs, intermediate and final outputs)](#per-module-reference)
   - [fastclean (PipeModules/FastClean.py)](#fastclean)
   - [kraken (PipeModules/Kraken.py)](#kraken)
   - [mapping (PipeModules/Mapping.py)](#mapping)
   - [coverage (PipeModules/Coverage.py)](#coverage)
   - [calling (PipeModules/Calling.py)](#calling)
   - [annotation_filter (PipeModules/AnnotationFilter.py)](#annotation_filter)
   - [consensus (PipeModules/Consensus.py)](#consensus)
   - [resistance (PipeModules/Resistance.py)](#resistance)
   - [typing (PipeModules/Typing.py)](#typing)
   - [distances (PipeModules/Distances.py)](#distances)
   - [getclusters / clusters (PipeModules/Clusters.py)](#getclusters)
   - [multiqc (PipeModules/Multiqc.py)](#multiqc)
   - [organize (PipeModules/Organize.py)](#organize)
   - [helpers: Repository, History, Version](#helpers)
- [Imported data files and where they are used (explicit list)](#important)
- [Troubleshooting and notes](#troubleshooting-and-notes)

<a id="about"></a>
About
-----
ThePipeline3 is a modular pipeline written in Python 3.7 intended for bacterial (MTBC) whole-genome sequencing processing: read cleaning, taxonomic filtering, mapping, variant calling, annotation, consensus generation, typing and resistance prediction.

<a id="short-dependencies"></a>
Short Dependencies
------------------
- Python 3.7
- Python packages: pandas, PyVCF (vcf), pairsnp (vendored in `data/libs/pairsnp-python/` but may be installed), and other stdlib modules.
- System/bioinformatics tools (expected under `Programs/` by default and configurable from `data/Paths/programs_path`): bwa, samtools, picard (MarkDuplicates jar), fastp, seqtk, pigz, kraken (and kraken-report/translate), bedtools (genomeCoverageBed), gatk (Mutect2), varscan (jar), minos (Singularity image), snpEff (jar), snp-sites, qualimap, pigz, snp-sites, others. See `Programs/` in the repository for bundled binaries and `data/Configs/software_versions.txt` for recorded versions.

<a id="layout-and-configuration-files-imported-by-the-pipeline"></a>
Layout and configuration files imported by the pipeline
-----------------------------------------------------
- `data/Paths/programs_path` — read by `PipeModules.Repository.Programs()`; contains program name mapping (format: name=/<path-to-program>) used by `Programs()` to build the `programs` dictionary consumed by modules.
- `data/Paths/data_path` — read by `PipeModules.Repository.Data()`; contains keys like `reference` and `krakenDB` used by default when callers don't provide them.
- `data/Configs/fastp_default.config` — default `fastp` parameter list used by `FastClean.ReadConfig()` when `--config default_config` is used.
- `data/Configs/multiqc.config` — used by the multiqc module to configure MultiQC.
- `data/Configs/software_versions.txt` — used by `PipeModules.Version.version()` and written into `<prefix>.history` files by `PipeModules.History.UpdateHistory()`; contains pinned program version info for reproducibility.
- `data/annotation_H37Rv.csv`, `data/H37Rv.annotation_new.tsv`, `data/catalogWHO2023_pipeline.csv`, `data/resistance_positions.csv`, `data/snp_phylo_fixed.tsv` — used by resistance, consensus, typing and annotation modules. Exact usage is documented in the module sections below.

<a id="how-the-driver-works-thepipeline3"></a>
How the driver works (`ThePipeline3`)
-----------------------------------
The main script `ThePipeline3` exposes subcommands (subparsers). Each subcommand corresponds to a module in `PipeModules/`. The driver:

- parses arguments with argparse (helper `percentage_float()` validates percentage inputs),
- maps subcommands to functions in `PipeModules` via `RunSubcommand(args)`,
- sets up environment variables where necessary (for example `LD_LIBRARY_PATH` when invoking bedtools via `Coverage`),
- prints a decorative banner when called with no args.

Top-level subcommands (and where they dispatch):

- fastclean -> PipeModules.FastClean.FastClean(args)
- kraken -> PipeModules.Kraken.Kraken(args)
- mapping -> PipeModules.Mapping.Mapping(args)
- coverage -> PipeModules.Coverage.CalcCoverage(args)
- calling -> PipeModules.Calling.Calling(args)
- annotation_filter -> PipeModules.AnnotationFilter.FilterSnps(args)
- consensus -> PipeModules.Consensus.Consensus(args)
- resistance -> PipeModules.Resistance.Resistance(args)
- typing -> PipeModules.Typing.Typing()
- distances -> PipeModules.Distances.Distances(args)
- organize -> PipeModules.Organize.Organize()
- multiqc -> PipeModules.Multiqc.Multiqc(args)
- getclusters -> PipeModules.Clusters.GetClusters(args)

<a id="per-module-reference"></a>
<a id="fastclean"></a>
FASTCLEAN (PipeModules/FastClean.py)
------------------------------------
**Purpose**
- Clean raw FASTQ reads using fastp, produce compressed cleaned FASTQ, a JSON and a cleanlog.

CLI parameters (as in `ThePipeline3` parser):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-f`, `--fastq` | One or two input FASTQ file paths. One file = single-end, two files = paired-end. | (required) | `nargs="*"` in parser; required for fastclean. |
| `-c`, `--config` | Path to a fastp config file (one parameter per line). | `default_config` | If `default_config`, reads `data/Configs/fastp_default.config`. |
| `-t`, `--threads` | Number of threads passed to fastp. | `1` | Up to 16 threads. |
| `-v`, `--verbose` | Print the fastp command before running. | (flag) | No value; present/absent. |
| `--phred64` | Treat input as Phred+64 (convert to Phred+33). | (flag) | Adds `--phred64` to fastp when set. |
| `-p`, `--prefix` | Output filename prefix used for produced files. | (required) | Required parameter used to name outputs. |

Required inputs
- One or two FASTQ files (gzipped or uncompressed depending on paths provided). The pipeline will not check file extension; it forwards the filenames to fastp.

Intermediate files created
- {prefix}.fastp.html — fastp HTML report (created, but removed immediately by the script; see "Files removed").
- {prefix}.fastp.json — fastp JSON report (kept)
- {prefix}.cleanlog — STDOUT/STDERR of fastp captured to this file
- {prefix}.clean.fastq.gz (single-end) or {prefix}.P1.clean.fastq.gz and {prefix}.P2.clean.fastq.gz (paired-end) — the cleaned and gzipped FASTQs

Files removed by the module
- {prefix}.fastp.html — removed in code because MultiQC can re-generate visualizations and the JSON is preserved.

Final outputs
- cleaned FASTQs (suffixes as above), `{prefix}.fastp.json`, `{prefix}.cleanlog`.

Notes
- The config file provides extra switches (like `--cut_by_quality3`, `--cut_window_size=10`, `--cut_mean_quality=20`, `--length_required=50`, `--correction`) — each parameter must be on its own (non-comment) line in `fastp_default.config` or a user-provided file.

Example
```bash
python3 ThePipeline3 fastclean -f sample_A_R1.fastq.gz sample_A_R2.fastq.gz -p sample_A -t 4 --config default_config
```
Outputs: `sample_A.P1.clean.fastq.gz`, `sample_A.P2.clean.fastq.gz`, `sample_A.fastp.json`, `sample_A.cleanlog`.

<a id="kraken"></a>
KRAKEN (PipeModules/Kraken.py)
--------------------------------
**Purpose**
- Classify reads with Kraken, translate taxonomy labels, extract reads that match a provided taxonomic string (by default `Mycobacterium tuberculosis`), and produce contamination reports.

CLI parameters (as in `ThePipeline3` parser):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-f`, `--fastq` | Input FASTQ file(s) to classify (single or paired). | (required for classify/filter) | `nargs="*"` — required when using `--classify` or `--filter`. |
| `-c`, `--compressed` | Tell Kraken input FASTQs are gzip compressed. | (flag) | Adds `--gzip-compressed` to kraken invocation. |
| `--paired` | Treat input as paired reads. | (flag) | Adjusts kraken invocation for paired data. |
| `-m`, `--matching` | String to match in kraken labels for read extraction. | `Mycobacterium tuberculosis` | Used by `MakeReadList` to identify reads to keep. |
| `-p`, `--prefix` | Output filename prefix. | (required) | Required; used to name outputs. |
| `--db` | Path to a Kraken database. | `False` (use data_path krakenDB) | If not provided, `data/Paths/data_path`'s `krakenDB` is used. |
| `-t`, `--threads` | Number of threads for Kraken. | `1` | Integer. |
| `--preload` | Preload the Kraken database into RAM. | (flag) | May require very large RAM. |
| `--classify` | Run Kraken classification to create `{prefix}.kraken`. | (flag) | Produces raw kraken output. |
| `--filter` | Extract reads matching `--matching` into filtered FASTQs. | (flag) | Produces `{prefix}.filtered.fastq` or paired filtered FASTQs. |
| `--noclean` | Keep Kraken intermediate files instead of deleting them. | (flag) | Prevents removal of `{prefix}.kraken`, `{prefix}.labels`, `{prefix}.filtered.readlist`. |
| `-r`, `--report` | Generate Kraken report and parse contaminants tables. | (flag) | Produces `{prefix}.kraken.report` and parsed tables. |

Required inputs
- FASTQ(s) for `--classify`/`--filter` options. If `--db` not provided, a `krakenDB` path must exist in `data/Paths/data_path`.

Intermediate files created
- {prefix}.kraken — kraken raw output (one line per read with classification). Created when `--classify`.
- {prefix}.labels — translated kraken labels (created by `kraken-translate`); used by `MakeReadList`.
- {prefix}.filtered.readlist — read ids of matching reads (created by `MakeReadList` using fgrep and cut) — used by `PickReads`.
- {prefix}.filtered.fastq / {prefix}.P1.filtered.fastq and {prefix}.P2.filtered.fastq — extracted reads matching taxonomy (compressed with pigz to .gz). Created by `PickReads`.
- {prefix}.kraken.report — generated by `kraken-report` and parsed to contaminants tables
- {prefix}.nfilter — small summary text created by `CountReads` with TOTAL and MTBC counts.

Files removed by the module
- By default (unless `--noclean`), the pipeline removes `{prefix}.labels`, `{prefix}.filtered.readlist`, and `{prefix}.kraken` after producing filtered fastqs and reports.

Final outputs
- `{prefix}.P1.filtered.fastq.gz` and `{prefix}.P2.filtered.fastq.gz` (paired mode) or `{prefix}.filtered.fastq.gz` (single-end) when `--filter` used.
- `{prefix}.kraken.report`, `{prefix}.genus.contaminants`, `{prefix}.species.contaminants` (if `--report` used).
- `{prefix}.nfilter` (read count summary)

Notes and pitfalls
- The `MakeReadList` implementation uses subprocess piping and decoding; ensure locale and encoding work when running on different systems.
- Kraken databases are very large — `--preload` will try to load the DB into RAM and may require 100s of GB.

Example
```bash
ThePipeline3 kraken -f sample_A.P1.clean.fastq.gz sample_A.P2.clean.fastq.gz --paired --compressed -p sample_A --classify --filter --report -t 5
```
Outputs: `sample_A.P1.filtered.fastq.gz`, `sample_A.P2.filtered.fastq.gz`, `sample_A.genus.contaminants`, `sample_A.species.contaminants`, `sample_A.nfilter`.

<a id="mapping"></a>
MAPPING (PipeModules/Mapping.py)
---------------------------------
**Purpose**
- Map reads to a reference with `bwa mem`, filter by mapQ, sort and (optionally) remove duplicates with Picard, hardclip reads with `samclip_h`, run QualiMap for QC and optionally convert to CRAM.

CLI parameters (as in `ThePipeline3`):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-f`, `--fastq` | One or two FASTQ files for mapping (single or paired). | (required) | `nargs="*"` in parser. |
| `-r`, `--reference` | Reference FASTA path to map reads against. | from `data/Paths/data_path` | If omitted, uses MTBC ancestor reference from data_path. |
| `-p`, `--prefix` | Output prefix for BAM/CRAM and QC files. | (required) | Used to name alignment outputs. |
| `-t`, `--threads` | Number of threads for bwa/samtools/qualimap. | `1` | Integer. |
| `--keep-dupbam` | Keep duplicate-marked BAMs alongside deduplicated BAMs. | (flag) | When set, original dup-marked files are preserved. |
| `-i`, `--index` | Run `bwa index -a bwtsw` on the provided reference. | (flag) | Use when first indexing a new reference. |
| `--no-dedup` | Skip duplicate removal (Picard MarkDuplicates). | (flag) | Avoids running dedup step. |
| `-mapq`, `--mapq` | MAPQ cutoff; filter alignments with MAPQ >= cutoff. | `0` | If >0, an `awk` filter is applied. |
| `-nhc`, `--no_hard_clipping` | Skip the hard-clipping step executed by `samclip_h`. | (flag) | Useful if hard clipping is not desired. |
| `-c`, `--cram` | Produce CRAM output (`{prefix}.sort.cram`) and remove BAM. | (flag) | Recommended to save space. Requires samtools with CRAM support and reference available. |

Required inputs
- FASTQ(s) and reference. If `--index` used, reference will be indexed by `bwa`.

Intermediate files created
- {prefix}.sort.bam — final sorted BAM after bwa mem and samtools sort (or deduped/renamed versions, see below)
- {prefix}.nd.sort.bam — intermediate output from Picard MarkDuplicates before renaming
- {prefix}.dup.metrix — metrics file produced by MarkDuplicates
- {prefix}_hardclipping.sort.bam — intermediate produced by hardclipping step (if executed)
- Qualimap output files (various) created in current dir
- {prefix}.sort.cram — (if `--cram` option chosen) CRAM produced with samtools view -C

Files removed / renamed by the module
- If duplicates are removed and `--keep-dupbam` is not set, `{}.nd.sort.bam` is renamed to `{}.sort.bam` (overwriting original). If `--keep-dupbam` is set, the original `.sort.bam` is saved as `{prefix}.dup.sort.bam` and the deduped file becomes `{prefix}.sort.bam`.
- If `--cram` chosen and CRAM conversion succeeds the script removes `{}.sort.bam`.
- The `hardclipping` function removes the original `{}.sort.bam` and replaces with the hardclipped renamed file.

Final outputs
- Alignment files: `{prefix}.sort.bam` (or `{prefix}.sort.cram` if `--cram`), `index` files as appropriate.
- Qualimap reports.

Notes
- Read group (@RG) is set using the prefix and some fixed fields to ensure GATK compatibility.

Example
```bash
ThePipeline3 mapping -f sample_A.P1.filtered.fastq.gz sample_A.P2.filtered.fastq.gz -p sample_A -t 8 -c
```
Outputs: `sample_A.sort.bam` (or `sample_A.sort.cram` if `-c`), `sample_A.dup.metrix`, QualiMap outputs.

<a id="coverage"></a>
COVERAGE (PipeModules/Coverage.py)
----------------------------------
**Purpose**
- Calculate per-base coverage with bedtools `genomeCoverageBed` and create summary files. Optionally filter samples that do not meet coverage thresholds and move their files to `NoPassCov/`.

CLI parameters (as in `ThePipeline3`):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-r`, `--reference` | Reference FASTA for genome length and bedtools. | data reference (from data_path) | If omitted, pipeline uses default data reference. |
| `-e`, `--extension` | Input alignment extension to use (`.sort.bam` or `.sort.cram`). | `bam` | Choices: `bam`, `cram`. |
| `-p`, `--prefix` | Sample prefix for coverage outputs. | (required) | Coverage computed for `{prefix}.sort.bam` or `{prefix}.sort.cram`. |
| `-f`, `--filter` | Move samples that don't meet thresholds to `NoPassCov/`. | (flag) | When set performs gating based on thresholds. |
| `--min-depth-mean` | Minimum mean depth required for pass. | `0` | Numeric. |
| `--min-depth-median` | Minimum median depth required for pass. | `20` | Numeric. |
| `--min-coverage` | Minimum fraction of genome covered at `--depth4cov`. | `0.95` | Validated via `percentage_float` (0.0-1.0). |
| `--keep-coverage` | Keep the per-base `{prefix}.coverage` files after summary. | (flag) | By default per-base coverage is removed. |
| `--depth4cov` | Depth threshold to consider a base as covered for genome coverage fraction. | `10` | Numeric depth threshold. |

Required inputs
- `{prefix}.sort.bam` or `{prefix}.sort.cram` depending on `--extension`. If `cram` selected, `CalcCoverage()` first uses samtools view to stream BAM to genomeCoverageBed.

Intermediate files created
- `{prefix}.coverage` — per-base coverage (created by genomeCoverageBed); huge file with one line per base (CHROM POS COV)
- `{prefix}.meancov` — summary file containing mean, median and genome_coverage computed using `depth4cov` threshold.

Files removed by the module
- `{prefix}.coverage` — removed unless `--keep-coverage` (because per-base files can be large).

Final outputs
- `{prefix}.meancov` — one-line human-readable summary, used by filters and reporting.

Example
```bash
ThePipeline3 coverage -p sample_A -e cram
```
Outputs: `sample_A.coverage` (unless `--keep-coverage` omitted), `sample_A.meancov`. If sample fails thresholds and `--filter` provided, files moved to `NoPassCov/`.

<a id="calling"></a>
CALLING (PipeModules/Calling.py)
--------------------------------
**Purpose**
- Variant discovery and adjudication: VarScan2 (pileup2snp/indel), GATK Mutect2-based calling, filtering of variants called in non-reliable genomic regions (PE/PPEs, highly repetitive regions), multiallelic handling, Minos adjudication (via Singularity) and SnpEff annotation (when reference corresponds to MTB ancestor/H37Rv). Variants with a frequence between 10%-90% are saved in a DR file (later used for resistance prediction) and fixed variants (freq>90%) are saved in a EPI file, which is further filtered by density.

CLI parameters (as parsed in `ThePipeline3`):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-e`, `--extension` | Extension used to find mapping files (e.g. `.cram` or `.sort.bam`). | `.sort.bam` | Validated by `check_extension()`. `.cram` files are coverted to `.sort.bam` files by default|
| `-r`, `--reference` | Reference FASTA path. | data reference (from data_path) | If omitted, pipeline uses default data reference, i.e. the MTB_ancestor reference, which originally comes from H37Rv but corrects false reference nucleotides in H37Rv. |
| `-p`, `--prefix` | Sample prefix. | (required) | Used to find mapping files and name outputs. |
| `-t`, `--threads` | Threads for GATK/samtools. | `1` | Integer. |
| `-kgvcf`, `--keep_gvcf` | Keep the generated `.gvcf` file. | (flag) | Otherwise the gVCF is removed by default. |
| `-kbam`, `--keep_bam` | Keep BAM files instead of converting to CRAM and removing them. | (flag) | Useful for debugging or workflow compatibility. |
| `-kmvcf`, `--keep_mutect2_vcf` | Keep Mutect2 intermediate VCF files. | (flag) | Controls whether Mutect2 VCFs are removed after their conversion to tabular format. |
| `-kna`, `--keep_not_annof` | Keep not-filtered-by-annotation VarScan and Mutect2 files | (flag) | |
| `-se`, `--single_end` | Indicate reads are single-end for Minos adjudication. | (flag) | Alters `Minos()` `--reads` invocation. |
| `-min_d`, `--min_depth` | Minimum depth considered callable. | `3` | Passed to Mutect2 `--callable-depth` and used in filters. |
| `-min_q`, `--min_qual` | Minimum base quality for callable positions. | `15` | Passed as Mutect2 `-mbq`. |
| `-min_f`, `--min_freq` | Minimum allele frequency to consider an alternative. | `0.05` | Validated via `percentage_float`. |
| `-filt_d`, `--filter_depth` | Depth threshold for EPI SNP filtering. | `20` | Numeric. |
| `-filt_f`, `--filter_freq` | Frequency threshold to include a SNP in EPI files. | `0.9` | Float between 0 and 1. |
| `--skip_dens_filt` | Skip density filtering stage. | (flag) | When set, densityfilter is not applied. |
| `-w`, `--window` | Sliding window size (bp) for density filtering. | `10` | Integer bp window. |
| `-dens`, `--density` | Density cutoff: number of SNPs within `window` that triggers removal. | `3` | Integer SNP count. |

Required inputs
- `{prefix}{ext}` where `ext` matches `--extension` (e.g. `.sort.bam` or `.sort.cram`) — mapping file produced by mapping module.
- `data/Paths/data_path` should provide default `reference` and other resources unless explicitly provided.
- Programs required (from `Programs()`): varscan (VarScan jar), gatk (gatk executable), minos (Singularity image) and snpEff (jar), samtools, genomeCoverageBed (bedtools), snp-sites etc.

High-level pipeline and generated intermediate files (ordered)

Additional intermediate files (created and removed during a standard `calling` run)
-----------------------------------------------------------------
The `calling` workflow now creates several short-lived intermediate files during processing. Most of these are removed in the final cleanup step (unless you pass flags like `--keep_gvcf`, `--keep_bam`, `--kmvcf`, or `--keep_not_annof`). For completeness, here is a compact list of common intermediate files you may see in a run:

- `{prefix}_mapq60.sort.bam` — MAPQ>=60 filtered BAM produced for VarScan (removed at end).
- `{prefix}.mpileup.remove` — samtools mpileup temporary file used by VarScan (removed after VarScan).
- `{prefix}.snp` — VarScan tabular output (renamed to `{prefix}.snp.varscan` during cleanup).
- `{prefix}.parsed.vcf` — VarScan-converted VCF for Minos (original saved as `.parsed.vcf.original_no_annoF` if annotation filtering is applied).
- `{prefix}_unfiltered.vcf`, `{prefix}.vcf`, `{prefix}_no_orientation.vcf` — Mutect2 intermediate VCFs produced during filtering and orientation removal (removed at end).
- `{prefix}.gvcf` and `{prefix}.gvcf.stats` — Mutect2 gVCF (kept only with `--keep_gvcf`).
- `{prefix}.indel.vcf`, `{prefix}.indel_sin_multiallelic.vcf`, `{prefix}.remade.indel.vcf` — intermediate indel extraction and filtering files (remade indels kept for downstream steps; originals may be removed).
- `{prefix}.coverage`, `{prefix}.lowcov.tsv` — bedtools per-base coverage and low-coverage positions (per-base coverage removed by default; lowcov TSV used to produce `{prefix}.lowcov`).
- `{prefix}.snp.vcf`, `{prefix}.multiallelic.snp.vcf`, `{prefix}.remade.snp.vcf` — Mutect2 SNP VCF and its multiallelic-splitted/remade versions. `{prefix}.remade.snp.vcf` is later renamed to `{prefix}.snp.mutect` in cleanup.
- `{prefix}_minos/` (directory) and `{prefix}.final_sin_wt.vcf` — Minos adjudication folder and extracted non-WT VCF (the folder is removed; final_sin_wt.vcf is removed during cleanup unless annotated output is kept).
- `{prefix}.final_sin_wt.vcf_complemented` — complemented Minos VCF with extra positions from VarScan/Mutect2 (temporary, removed after tabulation).
- `{prefix}.final_sin_wt_complemented_annotSnpEff.v` — SnpEff-annotated Minos VCF (created only when reference matches MTB_anc/H37Rv).
- `{prefix}.minos.raw.tab`, `{prefix}.filtered.minos.raw.tab` — tabular Minos+caller summary and its filtered version; the filtered file is renamed to `{prefix}.snp.minos` during cleanup.
- `{prefix}.snp.mutect.tab` — Mutect2 tabular conversion of `{prefix}.snp.mutect` (mutect VCF -> tab); original VCF removed unless `--kmvcf` set.
- `{prefix}.EPI.snp.nodensityfilter`, `{prefix}.dens_removed_snps.tab` — intermediate EPI/density filtering files (removed or consolidated into final `.EPI.snp.final.annoF`).

Note: the exact set of temporary files depends on the options you pass to `ThePipeline3 calling` (e.g., `--keep_gvcf`, `--kmvcf`, `--keep_not_annof`, `--keep_bam`). The full cleanup block in `PipeModules/Calling.py` shows the files that are explicitly removed or renamed at the end of the run.
1. **VarScan**
   - If `.cram` is input, converts to BAM first and sets `ext` to `.sort.bam`.
   - Creates `{prefix}_mapq60.sort.bam` (filtering MAPQ >= 60) using `mapq60_filter()`. This mapq60-filtered BAM will be the VarScan input. **In Mutect2, we keep the original because it already has its own internal filters.**
   - Runs `samtools mpileup` on the filtered BAM and writes `{prefix}.mpileup.remove`.
   - Runs VarScan `pileup2snp` producing `{prefix}.snp` (VarScan tab output).
   - Converts VarScan output to `{prefix}.parsed.vcf` (a small VCF-like file, created by the conversion code) so it can work as an input for Minos.
   - Removes `{prefix}.mpileup.remove` after VarScan.

2. **Mutect2**
   - Indexes BAM: `samtools index` creating `{prefix}.sort.bam.bai`.
   - Runs `gatk Mutect2` producing `{prefix}_unfiltered.vcf` (then filters orientation artifacts and produces `{prefix}.vcf`)
   - Produces gVCF `{prefix}.gvcf` using Mutect2 with `-ERC BP_RESOLUTION`.
   - SelectVariants splits SNPs/INDELs: `{prefix}.indel.vcf` and `{prefix}.snp.vcf`.
   - Extracts single-allele indels: `{prefix}.indel_sin_multiallelic.vcf` then filters to `{prefix}.remade.indel.vcf` (keeps indels with AF>=0.05 and DP>=10).
   - Calculates coverage with bedtools writing `{prefix}.coverage` and derives `{prefix}.lowcov.tsv` (positions with coverage < min_depth) and `{prefix}.lowcov` (converted later by `get_DR`).

3. **Annotation filtering & multiallelic handling**
   - `annoF_varscan_mutect()` uses `PipeModules.AnnotationFilter.LoadAnnotation()` loading `data/H37Rv.annotation_new.tsv` to mark positions to discard and writes `{prefix}.parsed.vcf_annoF`, `{prefix}.snp_annoF` and `{prefix}.snp.vcf_annoF` which are then renamed to replace the originals `{prefix}.parsed.vcf`, `{prefix].snp` and `{prefix}.snp.vcf`. It also renames originals to `.original_no_annoF`.
   - `filter_multiallelic_from_mutect2_snp()` reads `{prefix}.snp.vcf` and splits multiallelic records, saving multiallelic positions to `{prefix}.multiallelic.snp.vcf` and writing a remade VCF `{prefix}.remade.snp.vcf`, which contains multiallelic variants (when SNPs) but deployed in different rows so Minos can process the file.

4. **WT calling and Minos adjudication**
   - `callWT()` reads `{prefix}.gvcf` and writes `{prefix}.wt` (WT positions and depths), only saving positions with depth >= 3 and freq >= 90%

5. **Minos adjudication**
   - `Minos()` runs singularity-minos adjudication, producing a `{prefix}_minos` folder with Minos outputs including `final.vcf`; the code then extracts non-WT lines to `{prefix}.final_sin_wt.vcf` and removes the `{prefix}_minos` directory.
   - Annotates the Minos final VCF with SnpEff (produces `{prefix}.final_sin_wt.annotSnpEff.vcf`) if the reference is MTB_anc or H37Rv.

5. Tabulation and multi-caller reconciliation
   - `minos_raw_vcf_to_tab()` transforms Minos VCF to `{prefix}.minos.raw.tab`, complements Minos output with variants called by both VarScan and Mutect2 but not kept by Minos, computes mean frequencies and depths from VarScan and Mutect2, and writes `{prefix}.filtered.minos.raw.tab` after `filter_raw_minos()` which keeps only rows with depth >= 3 and freq >= 5%.

6. EPI, density filters and DR extraction
   - `filter_EPI()` produces `{prefix}.EPI.snp.nodensityfilter` and then (after densityfilter) produces `{prefix}.EPI.snp.final.annoF`.
   - `densityfilter()` writes `{prefix}.dens_removed_snps.tab` with the positions removed by density, and final `.EPI.snp.final.annoF` with the kept positions.
   - `get_DR()` reads `{prefix}.filtered.minos.raw.tab` to produce `{prefix}.DR.snp.final`, including gene annotations, and `{prefix}.lowcov`.

7. Cleanup and renames (the script's final steps)
   - Removes: `{prefix}.EPI.snp.nodensityfilter`, `{prefix}.final_sin_wt.vcf`, `{prefix}.lowcov.tsv`, `{prefix}.minos.raw.tab`, `{prefix}.parsed.vcf`, `{prefix}.snp.vcf`, `{prefix}.vcf`, `{prefix}.vcf.stats` (explicit removes in code). Also removes `{prefix}.snp.mutect` unless `-kmvcf` used.
   - Renames: `{prefix}.remade.snp.vcf` -> `{prefix}.snp.mutect`, `{prefix}.snp` -> `{prefix}.snp.varscan`, `{prefix}.filtered.minos.raw.tab` -> `{prefix}.snp.minos`.
   - Converts BAM -> CRAM and removes original `.sort.bam` unless `--keep_bam`.

Final outputs (user-facing)
- `{prefix}.DR.snp.final` — final annotated drug-relevant SNP list (tabular with Chrom, Position, Ref, Cons, VarFreq, Cov_total, VarAllele, Gene, Change, Ann, CodonWT, CodonVAR, Comments)
- `{prefix}.EPI.snp.final.annoF` — tabular SNPs used for population / epidemiology analyses
- `{prefix}.snp.minos`, `{prefix}.snp.varscan`, `{prefix}.snp.mutect(.tab)` (depending on flags) — per-caller outputs or reconciled outputs
- `{prefix}.wt` — wild-type callable positions
- `{prefix}.lowcov` — low coverage positions annotated with gene
- `{prefix}.gvcf` (kept if `--keep_gvcf`)

Example
```bash
ThePipeline3 calling -p sample_A -t 4 -e .cram
```
Major outputs: `sample_A.DR.snp.final`, `sample_A.EPI.snp.final.annoF`, `sample_A.snp.minos`, `sample_A.snp.varscan`, `sample_A.snp.mutect` (depending on flags), `sample_A.wt`, `sample_A.lowcov`, `{prefix}.gvcf` (if kept via `-kgvcf`).

Important flags and behaviour notes
- `--keep_gvcf` keeps the large `{prefix}.gvcf` file (otherwise removed) and toggles gzip compression when kept.
- `--keep_bam` keeps BAM files; by default pipeline converts to CRAM and removes BAM and its index.
- `--kmvcf` keeps the raw Mutect2 VCF; otherwise it is converted to tab and removed.

<a id="annotation_filter"></a>
ANNOTATION FILTER (PipeModules/AnnotationFilter.py)
-------------------------------------------------
**Purpose**
- Remove SNPs located in regions marked as DISCARD in `data/H37Rv.annotation_new.tsv` (these are typically phage, PE/PPE and repetitive regions). The filter produces filtered EPI SNP files and filtered VCFs.

CLI parameters
- Called by driver as `annotation_filter` and expects `-p/--prefix` argument (see ThePipeline3 parser). The function `FilterSnps(args)` uses `{prefix}.EPI.snp.final` and `{prefix}.EPI.snp.vcf` and writes `{prefix}.EPI.snp.final.annoF` and `{prefix}.EPI.snp.vcf.annoF`.

Inputs
- `{prefix}.EPI.snp.final`, `{prefix}.EPI.snp.vcf` (produced by Calling pipeline)
- `data/H37Rv.annotation_new.tsv` (used via `LoadAnnotation()`)

```bash
ThePipeline3 annotation_filter -s <input_file>
```

Outputs
- `{prefix}.EPI.snp.final.annoF` and `{prefix}.EPI.snp.vcf.annoF` — filtered versions; used downstream by consensus and distances modules.

<a id="consensus"></a>
CONSENSUS (PipeModules/Consensus.py)
-----------------------------------
**Purpose**
- Build non-redundant SNP tables across all samples, generate per-sample consensus FASTAs and concat them in a multiFASTA, run snp-sites to create SNP-only alignments, and produce final multi-FASTA and SNP tables useful for phylogenetics.

CLI parameters (driver `consensus` parser):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-i`, `--include` | Paths or folders containing sample `*.EPI.snp.final.annoF` files to include. | `.` | `nargs="*"` — one or more paths. |
| `-l`, `--sample_list` | Treat `paths` as a list of sample prefixes rather than folders. | (flag) | When set, `paths` are interpreted as sample prefixes. |
| `-t`, `--threads` | Number of threads for parallel steps. | `1` | Integer. |
| `-p`, `--prefix` | Output prefix for final multifasta/SNP table files. | (defaults to YYYYMMDD) | `dest`=`outfile` in parser. |

Inputs produced by the Calling pipeline:
- `{prefix}.EPI.snp.final.annoF` 
- `{prefix}.final_sin_wt.annotSnpEff.vcf` 
- `{prefix}.wt` 
- `{prefix}.indel.vcf`
- `{prefix}.snp.minos`
- `{prefix}.lowcov`
- `{prefix}.snp.varscan`
- `{prefix}.snp.mutect.tab`

Other inputs:
- `data/H37Rv.annotation_new.tsv` used for annotation context
- `data/resistance_positions.csv` used for annotating resistance positions in SNP table generation

Intermediate files and outputs
- `positions_total` and `snpeff_concat` — temporary concatenations used to build the SNP table (deleted later)
- `{outfile}.SNP_table.txt` — consolidated non-redundant SNP table with annotation and gene-level information (final)
- Per-sample `{prefix}.fas` consensus FASTA files (created per sample while generating individual FASTAs), then concatenated into `{outfile}.mf.fasta`, gapped version `{outfile}.mf_gap.fasta`, snp-sites FASTA `{outfile}.mf_gap.snp-sites.fasta` and a derived SNP table without variants associated with resistance `{outfile}.SNP_table.snp-sites.no-resis.txt`.

Cleanup
- Individual per-sample `.fas` files are removed after multi-FASTA generation.

Example
```bash
python3 ThePipeline3 consensus -i ./ -t 4 -p study_2025
```
Outputs: `study_2025.SNP_table.txt`, `study_2025.mf.fasta`, `study_2025.mf_gap.fasta`, `study_2025.mf_gap.snp-sites.fasta`, `study_2025.mf_gap.snp-sites.no-resis.fasta` and related SNP tables.

<a id="resistance"></a>
RESISTANCE (PipeModules/Resistance.py)
-------------------------------------
**Purpose**
- Generate per-sample resistance report files `{prefix}.res` from `{prefix}.DR.snp.final` using the WHO catalog `data/catalogWHO2023_pipeline.csv` and produce `resistance_report.csv` aggregating all samples.

CLI parameters
- -p, --prefix — if provided, analyze a single sample and produce `{prefix}.res`.
- -r, --report — if set, aggregate all `.res` files into `resistance_report.csv` in the working directory.

Inputs
- `{prefix}.DR.snp.final`, `{prefix}.lowcov`, `{prefix}.wt`, `{prefix}.remade.indel.vcf` — produced by Calling pipeline.
- `data/catalogWHO2023_pipeline.csv` — WHO mutations catalog (csv) used to map variants to drugs and confidence levels.

Example
```bash
ThePipeline3 resistance -p sample
ThePipeline3 resistance -r 
```

Outputs
- `{prefix}.res` — per-sample CSV-like file with Gene,Locus,Pos,CodonREF,CodonALT,Change,Freq,Drug,Confidence
- `resistance_report.csv` — aggregated report for all samples in folder when `CreateReport()` is called

Notes
- `CorrectRes()`, included in the script by default, performs corrections for triallelic/double/triple codon handling and rewrites a corrected `.res` file.

<a id="typing"></a>
TYPING (PipeModules/Typing.py)
------------------------------
**Purpose**
- Infer MTBC lineage from each `{prefix}.DR.snp.final` by matching phylogenetic marker positions in `data/snp_phylo_fixed.tsv` and produce `lineage_typing.csv` listing Sample, Infection type (clonal/mixed/undetermined) and typing.

Inputs
- `{prefix}.DR.snp.final` files found in working folder.
- `data/snp_phylo_fixed.tsv` — marker definitions.

Example
```bash
ThePipeline3 typing
```

Outputs
- `lineage_typing.csv` — lines of `Sample,Infection,Typing`.

<a id="distances"></a>
DISTANCES (PipeModules/Distances.py)
-----------------------------------
**Purpose**
- Calculate pairwise genetic distances.

CLI parameters (driver `distances` parser):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-p`, `--prefix` | Output prefix for distance files (`outfile`). | (defaults to YYYYMMDD) | `dest`=`outfile`. |
| `-f`, `--fasta` | Input FASTA file with sequences (e.g., `*.mf_gap.snp-sites.fasta`). | (required) | Required parameter. |
| `-t`, `--threads` | Number of threads for pairsnp. | `1` | Integer. |
| `-l`, `--limit` | Filter output to distances <= limit. | `-1` | If `-1` no filtering is applied. |

Example
```bash
ThePipeline3 distances -t 20 -f <multifasta> -p run1
```

Outputs
- `{outfile}.genetic_distances.tsv` — three-column table: Sequence_1, Distance, Sequence_2. If `--limit` set, an additional filtered file is written.

<a id="getclusters"></a>
GETCLUSTERS (PipeModules/Clusters.py)
-------------------------------------
**Purpose**
- Given a distances file (output of `distances`), compute clusters at a series of SNP thresholds (0,5,10,12,15 by default) and write `{prefix}.clusters_{threshold}.tsv` files (or `prov` intermediate) listing clusters and members.

CLI parameters (driver `getclusters` parser):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-d`, `--distances` | Input distances file (three-column from `distances`). | (required) | Required. |
| `-p`, `--prefix` | Output prefix for cluster files. | (required) | Required. |
| `-sep` | Separator for distances input (`Tab` or `Space`). | `Tab` | Controls parsing of the distances file. |
| `-osi` | Output clusters in alternate one-sample-per-line format. | (flag) | Produces OSI style output. What OSI means is an inside joke, don't worry about it, we just want to keep it that way |
| `-t`, `--threshold` | Compute clusters at a single threshold instead of the default set. | (none) | If omitted, computes default thresholds 0/5/10/12/15. |

Example
```bash
ThePipeline3 getclusters -d <distance_file> -thres 5 -p distance_file
```

Outputs
- `{outfile}.clusters_{threshold}.tsv` — cluster files filtered to clusters of size >=2.

<a id="multiqc"></a>
MULTIQC (PipeModules/Multiqc.py)
-------------------------------
**Purpose**
- Run MultiQC collecting logs and reports produced by the pipeline (fastp JSON, kraken reports, qualimap outputs, etc.). Uses a bundled MultiQC environment located in `Programs/MultiQC_env/` and config `data/Configs/multiqc.config`.

CLI parameters (driver `multiqc` parser):

| Parameter | Description | Default | Notes |
|---|---|---:|---|
| `-o`, `--output` | Output name for the MultiQC report. | `multiqc_report` | File/directory base name for the report. |
| `-f`, `--folder` | Folder to search for logs and reports to include. | `.` | Path to search recursively. |

Example
```bash
ThePipeline3 multiqc -o <output_name>
```

Outputs
- `multiqc` output directory and `multiqc_report.html` (or filename defined by `--output`).

<a id="organize"></a>
ORGANIZE (PipeModules/Organize.py)
---------------------------------
**Purpose**
- Move pipeline outputs into a set of user-facing folders (`reports`, `original_fastq`, `clean_fastq`, `filtered_fastq`, `BAM_CRAM`, `variants`, `resistance`) to tidy results.

Behaviour
- If a folder already exists, the script asks for confirmation interactively.
- Moves files according to filename patterns (see function body). It uses shell `mv` calls; if filenames conflict, files may be overwritten.

```bash
ThePipeline3 organize
```

<a id="helpers"></a>
Repository, History and Version helpers
-------------------------------------
**Repository**
- `PipeModules.Repository.Programs()` reads `{repo_root}/data/Paths/programs_path` and returns a dictionary mapping program names to absolute exe/jar paths. This is the central place to point the pipeline to local binaries.
- `PipeModules.Repository.Data()` reads `{repo_root}/data/Paths/data_path` and returns a dictionary with `reference`, `krakenDB` and other data resource paths used whenever modules don't receive explicit CLI paths.

History
- `PipeModules.History.UpdateHistory(command, program, prefix)` appends or creates `{cwd}/{prefix}.history` recording the program invoked, current working directory, timestamp, command string (or custom message) and version info from `data/Configs/software_versions.txt`.

Version
- `PipeModules.Version.version(program)` reads `data/Configs/software_versions.txt` and returns the recorded version string for a program key.

<a id="important"></a>
Imported data files and where they are used (explicit list)
--------------------------------------------------------
- `data/Paths/programs_path` — parsed by `PipeModules.Repository.Programs()` to find program paths.
- `data/Paths/data_path` — parsed by `PipeModules.Repository.Data()` to find default reference and krakenDB.
- `data/Configs/fastp_default.config` — read by `PipeModules.FastClean.ReadConfig()` for default fastp parameters.
- `data/Configs/multiqc.config` — used by `Multiqc` to configure MultiQC.
- `data/Configs/software_versions.txt` — read by `Version.version()` and written in `{prefix}.history`.
- `data/H37Rv.annotation_new.tsv` and `data/annotation_H37Rv.csv` — used by annotation and DR extraction (`AnnotationFilter`, `Calling`, `AnnotatePos`, `Resistance`, `Consensus`).
- `data/catalogWHO2023_pipeline.csv` — resistance catalog used by `Resistance.DetectResistance()` and `CorrectRes()`.
- `data/resistance_positions.csv` — used when building SNP_table/resistance tables.
- `data/snp_phylo_fixed.tsv` — phylogenetic SNP markers used by `Typing.Typing()`.

<a id="troubleshooting-and-notes"></a>
Troubleshooting and tips
------------------------
- Missing programs: ensure `data/Paths/programs_path` points to valid executables (or adjust PATH). Many heavy bioinformatics packages are expected in `Programs/`.
- Memory-heavy steps: Kraken DB classification and Minos adjudication (Singularity) can be resource intensive. Only run on machines with sufficient RAM/CPU.
- Encoding: Kraken `MakeReadList` uses `decode(sys.stdout.encoding)` — ensure a UTF-8 compatible locale to avoid decoding errors.
- Reproducibility: `{prefix}.history` files are written for each sample with the exact command that was run and version info from `data/Configs/software_versions.txt`.
