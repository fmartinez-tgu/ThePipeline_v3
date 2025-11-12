# ThePipeline3 — Quick start

This README gives a concise setup and example usage for the main modules
used in an analysis: `FastClean` (fastp), `Kraken` (classification),
`Coverage` (bedtools genomeCoverageBed) and `Calling` (VarScan2 / GATK
Mutect2 + Minos adjudication). The commands shown are dry-run examples for a
mock paired sample `sample_A` with input files:

- `sample_A_R1.fastq.gz` (read 1)
- `sample_A_R2.fastq.gz` (read 2)

Prerequisites
-------------

- Python 3.7 environment with `pandas`, `PyVCF` and `pairsnp` installed.
  See `environment.yml` (recommended) and `requirements.txt` (pip) added to
  the repository.
- Java (OpenJDK 8 or 11) available on PATH for GATK/VarScan/snpEff/picard.
- Command-line tools (install via conda/bioconda where possible):
  bwa, samtools, picard, fastp, kraken, seqtk, pigz, bedtools, varscan,
  snp-sites, snpEff, minos (containerized), snp-sites, qualimap, singularity.
- Many programs are bundled under `Programs/` in the repo. Ensure
  `data/Paths/programs_path` points to the correct relative paths or
  install tools system-wide.

Setup example (conda, recommended)
---------------------------------

Create the recommended environment and activate it:

```bash
conda env create -f environment.yml
conda activate thepipeline3
```

If you prefer pip, create a Python 3.7 venv and install python packages:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Step-by-step example (dry-run commands)
--------------------------------------

All commands assume you run them from the repository root and that
`ThePipeline3` is executable (`chmod +x ThePipeline3`). Replace `sample_A`
and fastq filenames as needed.

1) FastClean — clean reads with `fastp`

```bash
./ThePipeline3 fastclean -f sample_A_R1.fastq.gz sample_A_R2.fastq.gz -p sample_A -t 4
```

Outputs: `sample_A.P1.clean.fastq.gz`, `sample_A.P2.clean.fastq.gz`,
`sample_A.cleanlog`

2) Kraken — classify reads; optionally extract reads matching MTB

```bash
# classify, create report and filter reads matching default string (Mycobacterium tuberculosis)
./ThePipeline3 kraken -f sample_A.P1.clean.fastq.gz sample_A.P2.clean.fastq.gz --paired -p sample_A --classify --filter --report
```

Outputs: `sample_A.kraken`, `sample_A.labels`, `sample_A.filtered.readlist`,
`sample_A.P1.filtered.fastq.gz`, `sample_A.P2.filtered.fastq.gz`,
`sample_A.kraken.report`

3) Mapping — produce sorted BAM (required prior to Coverage and Calling)

```bash
./ThePipeline3 mapping -f sample_A.P1.filtered.fastq.gz sample_A.P2.filtered.fastq.gz -p sample_A -t 8
```

Output: `sample_A.sort.bam` (and index `sample_A.sort.bam.bai`)

4) Coverage — per-base coverage using bedtools

```bash
./ThePipeline3 coverage -e bam -p sample_A --min-depth-mean 0 --min-depth-median 20 --min-coverage 0.95 --depth4cov 10
```

Outputs: `sample_A.coverage` (per-base), `sample_A.meancov`

5) Calling — full variant calling pipeline (VarScan + Mutect2 + Minos etc.)

```bash
./ThePipeline3 calling -p sample_A -e .sort.bam -t 4 --min_d 3 --min_q 15 --min_f 0.05 --filt_d 20 --filt_f 0.9
```

Key outputs: `sample_A.snp.varscan`, `sample_A.snp.mutect.tab`,
`sample_A.remade.snp.vcf`, `sample_A.EPI.snp.final.annoF`, `sample_A.DR.snp.final`,
`sample_A.res`

Notes and tips
--------------
- If you want to use a specific reference, pass `-r /path/to/reference.fasta`
  to `mapping` and `calling` commands (otherwise the default defined in
  `data/Paths/data_path` is used).
- If `Programs/` includes binaries you prefer to use, ensure `data/Paths/programs_path`
  contains the correct (relative) paths; the pipeline's `Repository.Programs()`
  reads that file to construct tool paths.
- If you encounter failures during `calling` the history file `<prefix>.history`
  contains the exact commands executed and version information (see
  `PipeModules/History.py`). Use it to debug missing tools or wrong paths.
- Minos is executed via Singularity in this pipeline. Install Singularity on
  your host if you want to run Minos, or adapt the `Calling.Minos` call to use
  a local Minos installation.
  
 
