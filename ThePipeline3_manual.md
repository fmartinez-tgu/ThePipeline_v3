ThePipeline3 — Reference Manual
================================

This manual documents ThePipeline3 (Python-based), its submodules in `PipeModules/`, the command-line parameters exposed by `ThePipeline3` driver, the required inputs, intermediate files produced during execution, files removed by the pipeline, and the final outputs. It follows the style and structure of `Manual.md` but expands every module with parameter-by-parameter descriptions, examples and exact file flow so team members can run, debug and extend the pipeline with confidence.

Note: This document is generated from code in `ThePipeline3` and each `PipeModules/*.py` file. If you update the code, update this manual.

Contents
--------
- About
- Dependencies (short)
- Layout and configuration files imported by the pipeline
- How the driver works (`ThePipeline3`)
- Per-module reference (parameters, inputs, intermediate and final outputs)
  - fastclean (PipeModules/FastClean.py)
  - kraken (PipeModules/Kraken.py)
  - mapping (PipeModules/Mapping.py)
  - coverage (PipeModules/Coverage.py)
  - calling (PipeModules/Calling.py)
  - annotation_filter (PipeModules/AnnotationFilter.py)
  - consensus (PipeModules/Consensus.py)
  - resistance (PipeModules/Resistance.py)
  - typing (PipeModules/Typing.py)
  - distances (PipeModules/Distances.py)
  - getclusters / clusters (PipeModules/Clusters.py)
  - multiqc (PipeModules/Multiqc.py)
  - organize (PipeModules/Organize.py)
  - helpers: Repository, History, Version
- Examples (command lines)
- Troubleshooting and notes

About
-----
ThePipeline3 is a modular pipeline written in Python 3.7 intended for bacterial (MTBC) whole-genome sequencing processing: read cleaning, taxonomic filtering, mapping, variant calling, annotation, consensus generation, typing and resistance prediction.

Short Dependencies
------------------
- Python 3.7
- Python packages: pandas, PyVCF (vcf), pairsnp (vendored in `data/libs/pairsnp-python/` but may be installed), and other stdlib modules.
- System/bioinformatics tools (expected under `Programs/` by default and configurable from `data/Paths/programs_path`): bwa, samtools, picard (MarkDuplicates jar), fastp, seqtk, pigz, kraken (and kraken-report/translate), bedtools (genomeCoverageBed), gatk (Mutect2), varscan (jar), minos (Singularity image), snpEff (jar), snp-sites, qualimap, pigz, snp-sites, others. See `Programs/` in the repository for bundled binaries and `data/Configs/software_versions.txt` for recorded versions.

Layout and configuration files imported by the pipeline
-----------------------------------------------------
- `data/Paths/programs_path` — read by `PipeModules.Repository.Programs()`; contains program name mapping (format: name=/<path-to-program>) used by `Programs()` to build the `programs` dictionary consumed by modules.
- `data/Paths/data_path` — read by `PipeModules.Repository.Data()`; contains keys like `reference` and `krakenDB` used by default when callers don't provide them.
- `data/Configs/fastp_default.config` — default `fastp` parameter list used by `FastClean.ReadConfig()` when `--config default_config` is used.
- `data/Configs/multiqc.config` — used by the multiqc module to configure MultiQC.
- `data/Configs/software_versions.txt` — used by `PipeModules.Version.version()` and written into `<prefix>.history` files by `PipeModules.History.UpdateHistory()`; contains pinned program version info for reproducibility.
- `data/annotation_H37Rv.csv`, `data/H37Rv.annotation_new.tsv`, `data/catalogWHO2023_pipeline.csv`, `data/resistance_positions.csv`, `data/snp_phylo_fixed.tsv` — used by resistance, consensus, typing and annotation modules. Exact usage is documented in the module sections below.

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

Per-module reference (parameters, inputs, intermediate files, outputs)
------------------------------------------------------------------

Note on style: For every module I list: purpose, CLI parameters, required inputs, intermediate files created, files removed by the module, final outputs, and short comments about expected formats.

FASTCLEAN (PipeModules/FastClean.py)
------------------------------------
Purpose
- Clean raw FASTQ reads using fastp, produce compressed cleaned FASTQ, a JSON and a cleanlog.

CLI parameters (as in `ThePipeline3` parser):
- -f, --fastq (required, nargs="*") — one or two input fastq paths. If one file: single-end mode; if two: paired-end mode. Required.
- -c, --config (default: "default_config") — path to a config file listing additional fastp CLI parameters (one per line). If `default_config`, the pipeline reads `data/Configs/fastp_default.config`.
- -t, --threads (default: "1") — number of threads to pass to fastp. The code forces a maximum of 16 (fastp limitation).
- -v, --verbose (flag) — if set, prints the fastp command line before running it.
- --phred64 (flag) — if set, adds `--phred64` to fastp to indicate input is Phred+64; output will be Phred+33.
- -p, --prefix (required) — prefix used to name output files.

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

KRAKEN (PipeModules/Kraken.py)
--------------------------------
Purpose
- Classify reads with Kraken, translate taxonomy labels, extract reads that match a provided taxonomic string (by default `Mycobacterium tuberculosis`), and produce contamination reports.

CLI parameters (as in `ThePipeline3` parser):
- -f, --fastq (nargs="*") — input fastq(s) to classify. Required for `--classify` and `--filter` operations.
- -c, --compressed (flag) — whether input fastqs are gzip compressed; if set, `--gzip-compressed` will be added to kraken command.
- --paired (flag) — treat input as paired reads when running kraken.
- -m, --matching (default: "Mycobacterium tuberculosis") — string to match in labels for downstream filtering (used in `MakeReadList`).
- -p, --prefix (required) — prefix for output files.
- --db (default False) — path to a kraken database; if not provided, `PipeModules.Repository.Data()`'s `krakenDB` value is used.
- -t, --threads (default: "1") — number of threads for kraken.
- --preload (flag) — add `--preload` to kraken to pre-load the DB in RAM.
- --classify (flag) — run kraken classification and create `{prefix}.kraken`.
- --filter (flag) — extract reads matching the `--matching` pattern into `{prefix}.filtered.fastq` (.P1/.P2 for paired).
- --noclean (flag) — when present, do not delete intermediate kraken files (`{prefix}.kraken`, `{prefix}.labels`, `{prefix}.filtered.readlist`).
- -r, --report (flag) — generate `{prefix}.kraken.report` and parse to `{prefix}.genus.contaminants` and `{prefix}.species.contaminants`.

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

MAPPING (PipeModules/Mapping.py)
---------------------------------
Purpose
- Map reads to a reference with `bwa mem`, filter by mapQ, sort and (optionally) remove duplicates with Picard, hardclip reads with `samclip_h`, run QualiMap for QC and optionally convert to CRAM.

CLI parameters (as in `ThePipeline3`):
- -f, --fastq (required, nargs="*") — one or two fastq files (single or paired).
- -r, --reference (default: from `data/Paths/data_path`) — reference fasta path; when not provided, the pipeline uses the MTBC ancestor provided in data_path.
- -p, --prefix (required) — prefix for outputs.
- -t, --threads (default: "1") — number of threads for bwa/samtools/qualimap.
- --keep-dupbam (flag) — if set, the duplicate-marked BAMs are kept alongside deduplicated BAMs (renames are used by the code).
- -i, --index (flag) — run `bwa index -a bwtsw` on the reference before mapping (use only first time or when reference changed).
- --no-dedup (flag) — skip duplicate removal step.
- -mapq, --mapq (mapq_cutoff, default: 0) — integer threshold; if >0 an awk filter filters alignments to keep MAPQ >= mapq_cutoff.
- -nhc, --no_hard_clipping (flag) — if set, skip the `hardclipping` step executed by `samclip_h`.
- -c, --cram (flag) — if set, produce CRAM `.{prefix}.sort.cram` and remove `.sort.bam`.

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

COVERAGE (PipeModules/Coverage.py)
----------------------------------
Purpose
- Calculate per-base coverage with bedtools `genomeCoverageBed` and create summary files. Optionally filter samples that do not meet coverage thresholds and move their files to `NoPassCov/`.

CLI parameters (as in `ThePipeline3`):
- -r, --reference (default: data reference) — reference for genome length and bedtools.
- -e, --extension (default: "bam", choices {"bam","cram"}) — whether to read `.sort.bam` or `.sort.cram`.
- -p, --prefix (required) — sample prefix (the coverage is computed for `{prefix}.sort.bam` or `{prefix}.sort.cram`).
- -f, --filter (flag) — If set, samples that don't meet thresholds are moved to `NoPassCov/`.
- --min-depth-mean (default: 0) — minimum mean depth required (MeanCoverage uses this when filtering).
- --min-depth-median (default: 20) — minimum median depth required for pass.
- --min-coverage (default: 0.95) — minimum genomic fraction covered at `--depth4cov` (validated via `percentage_float` to be between 0.0 and 1.0).
- --keep-coverage (flag) — if set, keep `{prefix}.coverage` files; otherwise they are removed after meancov calculation.
- --depth4cov (default: 10) — minimum depth considered for counting a base as covered for genome coverage.

Required inputs
- `{prefix}.sort.bam` or `{prefix}.sort.cram` depending on `--extension`. If `cram` selected, `CalcCoverage()` first uses samtools view to stream BAM to genomeCoverageBed.

Intermediate files created
- `{prefix}.coverage` — per-base coverage (created by genomeCoverageBed); huge file with one line per base (CHROM POS COV)
- `{prefix}.meancov` — summary file containing mean, median and genome_coverage computed using `depth4cov` threshold.

Files removed by the module
- `{prefix}.coverage` — removed unless `--keep-coverage` (because per-base files can be large).

Final outputs
- `{prefix}.meancov` — one-line human-readable summary, used by filters and reporting.

CALLING (PipeModules/Calling.py)
--------------------------------
Purpose
- Variant discovery and adjudication: VarScan2 (pileup2snp/indel), GATK Mutect2-based calling, multiallelic handling, Minos adjudication (via Singularity), SnpEff annotation (when reference corresponds to MTB ancestor/H37Rv), filtering, EPI and density filtering, DR extraction and resistance-ready outputs.

CLI parameters (as parsed in `ThePipeline3`):
- -e, --extension (dest: ext, default: ".sort.bam") — extension used to find mapping files; should be `.sort.bam` or `.sort.cram` (Mutect2/VarScan expect either). The code validates via `check_extension()`.
- -r, --reference (default: data reference) — reference fasta path.
- -p, --prefix (required) — sample prefix.
- -t, --threads (default: "1") — threads for GATK/samtools.
- -kgvcf, --keep_gvcf (flag) — keep the generated `.gvcf` file (otherwise removed by default).
- -kbam, --keep_bam (flag) — keep BAM files instead of converting to CRAM and removing BAMs.
- -kmvcf, --keep_mutect2_vcf (flag) — keep the original Mutect2 VCF (mutect `*_unfiltered.vcf` / intermediate VCFs). Controls whether `mutect2_vcf_to_tab()` removes original VCFs.
- -se, --single_end (flag) — indicate the reads are single-end for Minos adjudication (changes Minos `--reads` invocation and file names).
- -min_d, --min_depth (default: "3") — minimum depth considered callable (passed to Mutect2 `--callable-depth` and used in other parts). String in parser but used as numeric later.
- -min_q, --min_qual (default: "15") — minimum base quality for callable positions (Mutect2 `-mbq` parameter).
- -min_f, --min_freq (default: "0.05") — minimum allele frequency to consider an alternative (applies to some filters/VarScan settings) — validated by `percentage_float`.
- -filt_d, --filter_depth (default: 20) — depth threshold used in EPI SNP filtering (EPI.snp file creation).
- -filt_f, --filter_freq (default: 0.9) — frequency threshold (as float 0.9) to include a SNP in EPI.snp files.
- --skip_dens_filt (flag) — skip density filter stage (densityfilter).
- -w, --window (default: 10) — sliding window size (bp) used for density filtering.
- -dens, --density (default: 3) — density cutoff: number of SNPs within `window` triggers removal.

Required inputs
- `{prefix}{ext}` where `ext` matches `--extension` (e.g. `.sort.bam` or `.sort.cram`) — mapping file produced by mapping module.
- `data/Paths/data_path` should provide default `reference` and other resources unless explicitly provided.
- Programs required (from `Programs()`): varscan (VarScan jar), gatk (gatk executable), minos (Singularity image) and snpEff (jar), samtools, genomeCoverageBed (bedtools), snp-sites etc.

High-level pipeline and generated intermediate files (ordered)
1. VarScan path
   - Creates `{prefix}_mapq60.sort.bam` (filtering MAPQ >= 60) using `mapq60_filter()`.
   - Runs `samtools mpileup` on the filtered BAM and writes `{prefix}.mpileup.remove`.
   - Runs VarScan `pileup2snp` producing `{prefix}.snp` (VarScan tab output).
   - Converts VarScan output to `{prefix}.parsed.vcf` (a small VCF-like file, created by the conversion code).
   - Removes `{prefix}.mpileup.remove` after VarScan.

2. Mutect2 path
   - If `.cram` is input, converts to BAM first and sets `ext` to `.sort.bam`.
   - Indexes BAM: `samtools index` creating `{prefix}.sort.bam.bai`.
   - Runs `gatk Mutect2` producing `{prefix}_unfiltered.vcf` (then filters orientation artifacts and produces `{prefix}.vcf`)
   - Produces gVCF `{prefix}.gvcf` using Mutect2 with `-ERC BP_RESOLUTION`.
   - SelectVariants splits SNPs/INDELs: `{prefix}.indel.vcf` and `{prefix}.snp.vcf`.
   - Extracts single-allele indels: `{prefix}.indel_sin_multiallelic.vcf` then filters to `{prefix}.remade.indel.vcf` (keeps indels with AF>=0.05 and DP>=10).
   - Calculates coverage with bedtools writing `{prefix}.coverage` and derives `{prefix}.lowcov.tsv` (positions with coverage < min_depth) and `{prefix}.lowcov` (converted later by `get_DR`).

3. Annotation filtering & multiallelic handling
   - `annoF_varscan_mutect()` uses `PipeModules.AnnotationFilter.LoadAnnotation()` loading `data/H37Rv.annotation_new.tsv` to mark positions to discard and writes `{prefix}.parsed.vcf_annoF`, `{prefix}.snp_annoF` which are then renamed to replace the originals. It also renames originals to `.original_no_annoF`.
   - `filter_multiallelic_from_mutect2_snp()` reads `{prefix}.snp.vcf` and splits multiallelic records, saving multiallelic positions to `{prefix}.multiallelic.snp.vcf` and writing a remade VCF `{prefix}.remade.snp.vcf`.

4. WT calling and Minos adjudication
   - `callWT()` reads `{prefix}.gvcf` and writes `{prefix}.wt` (WT positions and depths).
   - `Minos()` runs singularity-minos adjudication, producing a `{prefix}_minos` folder with Minos outputs including `final.vcf`; the code then extracts non-WT lines to `{prefix}.final_sin_wt.vcf` and removes the `{prefix}_minos` directory.
   - Annotates the Minos final VCF with SnpEff (produces `{prefix}.final_sin_wt.annotSnpEff.vcf`) if the reference is MTB_anc or H37Rv.

5. Tabulation and multi-caller reconciliation
   - `minos_raw_vcf_to_tab()` transforms Minos VCF to `{prefix}.minos.raw.tab`, computes mean frequencies and depths from VarScan and Mutect2, and writes `{prefix}.filtered.minos.raw.tab` after `filter_raw_minos()` which keeps only rows with depth >= 3 and freq >= 5%.

6. EPI, density filters and DR extraction
   - `filter_EPI()` produces `{prefix}.EPI.snp.nodensityfilter` and then (after densityfilter) produces `{prefix}.EPI.snp.final.annoF`.
   - `densityfilter()` writes `{prefix}.dens_removed_snps.tab` with the positions removed by density, and final `.EPI.snp.final.annoF` with the kept positions.
   - `get_DR()` reads `{prefix}.filtered.minos.raw.tab` and annotation (annotation CSV) to produce `{prefix}.DR.snp.final` and `{prefix}.lowcov` and other outputs.

7. Cleanup and renames (the script's final steps)
   - Removes: `{prefix}.EPI.snp.nodensityfilter`, `{prefix}.final_sin_wt.vcf`, `{prefix}.lowcov.tsv`, `{prefix}.minos.raw.tab`, `{prefix}.parsed.vcf`, `{prefix}.snp.vcf`, `{prefix}.vcf`, `{prefix}.vcf.stats` (explicit removes in code). Also removes `{prefix}.snp.mutect` unless `-kmvcf` used.
   - Renames: `{prefix}.remade.snp.vcf` -> `{prefix}.snp.mutect`, `{prefix}.snp` -> `{prefix}.snp.varscan`, `{prefix}.filtered.minos.raw.tab` -> `{prefix}.snp.minos`.
   - Converts BAM -> CRAM and removes original `.sort.bam` unless `--keep_bam`.
   - Removes `{prefix}_minos` directory after extracting data.

Final outputs (user-facing)
- `{prefix}.DR.snp.final` — final annotated drug-relevant SNP list (tabular with Chrom, Position, Ref, Cons, VarFreq, Cov_total, VarAllele, Gene, Change, Ann, CodonWT, CodonVAR, Comments)
- `{prefix}.EPI.snp.final.annoF` — tabular SNPs used for population / epidemiology analyses
- `{prefix}.snp.minos`, `{prefix}.snp.varscan`, `{prefix}.snp.mutect` (depending on flags) — per-caller outputs or reconciled outputs
- `{prefix}.wt` — wild-type callable positions
- `{prefix}.lowcov` — low coverage positions annotated with gene
- `{prefix}.gvcf` (kept if `--keep_gvcf`)

Important flags and behaviour notes
- `--keep_gvcf` keeps the large `{prefix}.gvcf` file (otherwise removed) and toggles gzip compression when kept.
- `--keep_bam` keeps BAM files; by default pipeline converts to CRAM and removes BAM and its index.
- `--kmvcf` keeps the raw Mutect2 VCF; otherwise it is converted to tab and removed.

ANNOTATION FILTER (PipeModules/AnnotationFilter.py)
-------------------------------------------------
Purpose
- Remove SNPs located in regions marked as DISCARD in `data/H37Rv.annotation_new.tsv` (these are typically phage, PE/PPE and repetitive regions). The filter produces filtered EPI SNP files and filtered VCFs.

CLI parameters
- Called by driver as `annotation_filter` and expects `-p/--prefix` argument (see ThePipeline3 parser). The function `FilterSnps(args)` uses `{prefix}.EPI.snp.final` and `{prefix}.EPI.snp.vcf` and writes `{prefix}.EPI.snp.final.annoF` and `{prefix}.EPI.snp.vcf.annoF`.

Inputs
- `{prefix}.EPI.snp.final`, `{prefix}.EPI.snp.vcf` (produced by Calling pipeline)
- `data/H37Rv.annotation_new.tsv` (used via `LoadAnnotation()`)

Outputs
- `{prefix}.EPI.snp.final.annoF` and `{prefix}.EPI.snp.vcf.annoF` — filtered versions; used downstream by consensus and distances modules.

CONSENSUS (PipeModules/Consensus.py)
-----------------------------------
Purpose
- Build non-redundant SNP tables across all samples, generate per-sample consensus FASTAs and multiFASTAs, run snp-sites to create SNP-only alignments, and produce final multi-FASTA and SNP tables useful for phylogenetics.

CLI parameters (driver `consensus` parser):
- -i, --include (dest `paths`, default '.') nargs="*" — one or more paths (folders) to samples that contain `*.EPI.snp.final.annoF` files; these will be included in the multifasta.
- -l, --sample_list (flag) — if set, treat `paths` as a list of prefixes (i.e., a list of sample prefixes) to include instead of paths with folders.
- -t, --threads (default 1) — number of threads for parallel steps.
- -p, --prefix (dest `outfile`) — prefix used to name final files (defaults to YYYYMMDD if not provided).

Inputs
- `{prefix}.EPI.snp.final.annoF` — annotated SNP files for each sample (present in each sample folder or prefix list)
- `{prefix}.final_sin_wt.annotSnpEff.vcf` — SnpEff annotated Minos outputs (concatenated and deduped for annotation information)
- `data/H37Rv.annotation_new.tsv` used for annotation context
- `data/resistance_positions.csv` used for annotating resistance positions in SNP table generation

Intermediate files and outputs
- `positions_total` and `snpeff_concat` — temporary concatenations used to build the SNP table (deleted later)
- `{outfile}.SNP_table.txt` — consolidated non-redundant SNP table with annotation and gene-level information (final)
- Per-sample `{prefix}.fas` consensus FASTA files (created per sample while generating individual FASTAs), then concatenated into `{outfile}.mf.fasta`, gapped version `{outfile}.mf_gap.fasta`, snp-sites FASTA `{outfile}.mf_gap.snp-sites.fasta` and derived SNP table files such as `{outfile}.SNP_table.snp-sites.no-resis.txt`.

Cleanup
- Individual per-sample `.fas` files are removed after multi-FASTA generation.

RESISTANCE (PipeModules/Resistance.py)
-------------------------------------
Purpose
- Generate per-sample resistance report files `{prefix}.res` from `{prefix}.DR.snp.final` using the WHO catalog `data/catalogWHO2023_pipeline.csv` and produce `resistance_report.csv` aggregating all samples.

CLI parameters
- -p, --prefix — if provided, analyze a single sample and produce `{prefix}.res`.
- -r, --report — if set, aggregate all `.res` files into `resistance_report.csv` in the working directory.

Inputs
- `{prefix}.DR.snp.final`, `{prefix}.lowcov`, `{prefix}.wt`, `{prefix}.remade.indel.vcf` — produced by Calling pipeline.
- `data/catalogWHO2023_pipeline.csv` — WHO mutations catalog (csv) used to map variants to drugs and confidence levels.

Outputs
- `{prefix}.res` — per-sample CSV-like file with Gene,Locus,Pos,CodonREF,CodonALT,Change,Freq,Drug,Confidence
- `resistance_report.csv` — aggregated report for all samples in folder when `CreateReport()` is called

Notes
- `CorrectRes()` performs corrections for triallelic/double/triple codon handling and rewrites a corrected `.res` file.

TYPING (PipeModules/Typing.py)
------------------------------
Purpose
- Infer MTBC lineage from each `{prefix}.DR.snp.final` by matching phylogenetic marker positions in `data/snp_phylo_fixed.tsv` and produce `lineage_typing.csv` listing Sample, Infection type (clonal/mixed/undetermined) and typing.

Inputs
- `{prefix}.DR.snp.final` files found in working folder.
- `data/snp_phylo_fixed.tsv` — marker definitions.

Outputs
- `lineage_typing.csv` — lines of `Sample,Infection,Typing`.

DISTANCES (PipeModules/Distances.py)
-----------------------------------
Purpose
- Calculate pairwise genetic distances using the vendored/installed `pairsnp` package.

CLI parameters (driver `distances` parser):
- -p, --prefix (dest `outfile`) — output prefix (defaults to YYYYMMDD).
- -f, --fasta (required) — FASTA file containing sequences (e.g., `*.mf_gap.snp-sites.fasta`).
- -t, --threads (default 1) — threads for pairsnp.
- -l, --limit (default -1) — optionally filter output distances to those <= limit.

Outputs
- `{outfile}.genetic_distances.tsv` — three-column table: Sequence_1, Distance, Sequence_2. If `--limit` set, an additional filtered file is written.

GETCLUSTERS (PipeModules/Clusters.py)
-------------------------------------
Purpose
- Given a distances file (output of `distances`), compute clusters at a series of SNP thresholds (0,5,10,12,15 by default) and write `{prefix}.clusters_{threshold}.tsv` files (or `prov` intermediate) listing clusters and members.

CLI parameters (driver `getclusters` parser):
- -d, --distances (required) — distances file.
- -p, --prefix (required) — output prefix.
- -sep (default Tab) — separator of distances input (Tab or Space).
- -osi (flag) — if set, produce alternate per-line sample output format (one sample per line) compatible with the 'Senyorito Irving' style.
- -t, --threshold — single threshold to compute (otherwise computes default set 0/5/10/12/15).

Outputs
- `{outfile}.clusters_{threshold}.tsv` — cluster files filtered to clusters of size >=2.

MULTIQC (PipeModules/Multiqc.py)
-------------------------------
Purpose
- Run MultiQC collecting logs and reports produced by the pipeline (fastp JSON, kraken reports, qualimap outputs, etc.). Uses a bundled MultiQC environment located in `Programs/MultiQC_env/` and config `data/Configs/multiqc.config`.

CLI parameters (driver `multiqc` parser):
- -o, --output (default `multiqc_report`) — output filename.
- -f, --folder (default '.') — folder to search for reports and logs.

Outputs
- `multiqc` output directory and `multiqc_report.html` (or filename defined by `--output`).

ORGANIZE (PipeModules/Organize.py)
---------------------------------
Purpose
- Move pipeline outputs into a set of user-facing folders (`reports`, `original_fastq`, `clean_fastq`, `filtered_fastq`, `BAM_CRAM`, `variants`, `resistance`) to tidy results.

Behaviour
- If a folder already exists, the script asks for confirmation interactively.
- Moves files according to filename patterns (see function body). It uses shell `mv` calls; if filenames conflict, files may be overwritten.

Repository, History and Version helpers
-------------------------------------
Repository
- `PipeModules.Repository.Programs()` reads `{repo_root}/data/Paths/programs_path` and returns a dictionary mapping program names to absolute exe/jar paths. This is the central place to point the pipeline to local binaries.
- `PipeModules.Repository.Data()` reads `{repo_root}/data/Paths/data_path` and returns a dictionary with `reference`, `krakenDB` and other data resource paths used whenever modules don't receive explicit CLI paths.

History
- `PipeModules.History.UpdateHistory(command, program, prefix)` appends or creates `{cwd}/{prefix}.history` recording the program invoked, current working directory, timestamp, command string (or custom message) and version info from `data/Configs/software_versions.txt`.

Version
- `PipeModules.Version.version(program)` reads `data/Configs/software_versions.txt` and returns the recorded version string for a program key.

Examples (quick start)
----------------------
These example commands are dry-run style CLI invocations. Replace `sample_A_R1.fastq.gz` and `sample_A_R2.fastq.gz` with your files and `PREFIX` with the sample name.

1) Clean reads with fastp (paired):

```bash
python3 ThePipeline3 fastclean -f sample_A_R1.fastq.gz sample_A_R2.fastq.gz -p sample_A -t 4 --config default_config
```

Outputs: `sample_A.P1.clean.fastq.gz`, `sample_A.P2.clean.fastq.gz`, `sample_A.fastp.json`, `sample_A.cleanlog`.

2) Classify reads and extract MTBC reads with Kraken (paired):

```bash
python3 ThePipeline3 kraken -f sample_A.P1.clean.fastq.gz sample_A.P2.clean.fastq.gz --paired -p sample_A --classify --filter --report
```

Outputs: `sample_A.P1.filtered.fastq.gz`, `sample_A.P2.filtered.fastq.gz`, `sample_A.genus.contaminants`, `sample_A.species.contaminants`, `sample_A.nfilter`.

3) Map filtered reads to reference:

```bash
python3 ThePipeline3 mapping -f sample_A.P1.filtered.fastq.gz sample_A.P2.filtered.fastq.gz -p sample_A -t 8 -c
```

Outputs: `sample_A.sort.bam` (or `sample_A.sort.cram` if `-c`), `sample_A.dup.metrix`, QualiMap outputs.

4) Calculate coverage and optionally apply coverage filter:

```bash
python3 ThePipeline3 coverage -p sample_A -e bam --min-depth-mean 10 --min-depth-median 30 --min-coverage 0.95
```

Outputs: `sample_A.coverage` (unless `--keep-coverage` omitted), `sample_A.meancov`. If sample fails thresholds and `--filter` provided, files moved to `NoPassCov/`.

5) Variant calling pipeline (VarScan/Mutect2/Minos) :

```bash
python3 ThePipeline3 calling -p sample_A -t 4 -min_d 3 -min_q 15 -min_f 0.05
```

Major outputs: `sample_A.DR.snp.final`, `sample_A.EPI.snp.final.annoF`, `sample_A.snp.minos`, `sample_A.snp.varscan`, `sample_A.snp.mutect` (depending on flags), `sample_A.wt`, `sample_A.lowcov`, `{prefix}.gvcf` (if kept via `-kgvcf`).

6) Generate consensus and multifasta

```bash
python3 ThePipeline3 consensus -i ./ -t 4 -p study_2025
```

Outputs: `study_2025.SNP_table.txt`, `study_2025.mf.fasta`, `study_2025.mf_gap.fasta`, `study_2025.mf_gap.snp-sites.fasta`, and related SNP tables.

Files the pipeline often removes (summary)
----------------------------------------
- FastClean: removes `{prefix}.fastp.html` (keeps JSON and cleanlog).
- Kraken: removes `{prefix}.kraken`, `{prefix}.labels`, `{prefix}.filtered.readlist` unless `--noclean`.
- Mapping: may remove `.sort.bam` when converting to `.sort.cram` (unless `--keep_bam`), and may rename files when deduping.
- Coverage: removes `{prefix}.coverage` unless `--keep-coverage`.
- Calling: removes many intermediate files at end; the code explicitly `os.remove()`s or renames a variety of files (see list in Calling() near the end). Important removed files include: `{prefix}.EPI.snp.nodensityfilter`, `{prefix}.final_sin_wt.vcf`, `{prefix}.lowcov.tsv`, `{prefix}.minos.raw.tab`, `{prefix}.parsed.vcf`, `{prefix}.snp.vcf`, `{prefix}.vcf`, `{prefix}.vcf.stats`, and possibly `{prefix}.snp.mutect` (depending on `-kmvcf` flag). The pipeline also deletes `{prefix}_minos` directory at the end.

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

Troubleshooting and tips
------------------------
- Missing programs: ensure `data/Paths/programs_path` points to valid executables (or adjust PATH). Many heavy bioinformatics packages are expected in `Programs/`.
- Memory-heavy steps: Kraken DB classification and Minos adjudication (Singularity) can be resource intensive. Only run on machines with sufficient RAM/CPU.
- Encoding: Kraken `MakeReadList` uses `decode(sys.stdout.encoding)` — ensure a UTF-8 compatible locale to avoid decoding errors.
- Reproducibility: `{prefix}.history` files are written for each sample with the exact command that was run and version info from `data/Configs/software_versions.txt`.

Next steps and maintenance notes
-------------------------------
- Keep `data/Configs/software_versions.txt` up-to-date when upgrading tools; it is also used for `history` tracking.
- Consider adding a `--dry-run` or `--print-cmds` wrapper (the code already prints fastp command with `--verbose`) that echoes each external command before running it — useful for audit.
- Add small unit tests or smoke datasets (tiny FASTQ) to let users validate environment installs.

Acknowledgements
----------------
Manual generated by programmatic inspection of `ThePipeline3` and `PipeModules/` files in this repository. Use this as the authoritative reference for current code; update whenever code changes.

End of ThePipeline3_manual.md
