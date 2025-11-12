# ThePipeline v.2.0 (ThePipeline2)

- [ThePipeline v.2.0 (ThePipeline2)](#thepipeline-v20-thepipeline2)
    - [About](#about)
    - [Dependencies](#dependencies)
    - [Structure](#structure)
      - [data](#data)
      - [Programs](#programs)
      - [PipeModules](#pipemodules)
    - [Instalation](#instalation)
    - [Updating ThePipeline2](#updating-thepipeline2)
- [Running ThePipeline2](#running-thepipeline2)
- [Modules](#modules)
  - [FASTQ Filtering (fastclean)](#fastq-filtering-fastclean)
  - [Taxonomic filtering (kraken)](#taxonomic-filtering-kraken)
  - [Mapping (mapping)](#mapping-mapping)
  - [Variant calling (calling)](#variant-calling-calling)
  - [Coverage calculation and filtering (coverage)](#coverage-calculation-and-filtering-coverage)
  - [Filter by annotation (annotation\_filter)](#filter-by-annotation-annotation_filter)
  - [Generate multifastas (consensus)](#generate-multifastas-consensus)
  - [Organizing results (organize)](#organizing-results-organize)
  - [Detect genomic determinants of resistance (resistance)](#detect-genomic-determinants-of-resistance-resistance)
  - [Lineage typing (typing)](#lineage-typing-typing)
  - [Calculate genetic distances (distances)](#calculate-genetic-distances-distances)
  - [Calculate transmission clusters (getclusters)](#calculate-transmission-clusters-getclusters)
  - [Generate MultiQC report (multiqc)](#generate-multiqc-report-multiqc)

### About

_ThePipeline_ is a bioinformatic analysis pipeline developed by members of the Tuberculosis Genomics Unit (TGU) at IBV-CSIC. It was initially created for analyzing Whole Genome Sequencing (WGS) data from _Mycobacterium tuberculosis_ complex (MTBC) samples. Nevertheless, it could potentially work on a wide range of bacterial species, after tuning the default parameters and options. _ThePipeline2_ is the new version of _ThePipeline_. Some functions remain without changes, some others have changed a bit, and some other are new. So please, read this manual to understand how each of the modules work and the different outputs and functionalities

_ThePipeline2_ has several modules for performing common bioinformatic analysis. Basic runs consist on the execution of the following sequential steps:

- Fastq filtering
- Taxonomic filtering
- Mapping to a reference
- Variant calling discovery and filtering

Each of these steps is executed by an independent module, which uses as an input the output of the previous one. In addition to these 'core' modules, _ThePipeline2_ has other modules with extended functions. See bellow for a complete description of all the modules, including its main arguments and options, and the output obtained.


### Dependencies
_ThePipeline_ depend on several programs that can be found in the <u>ThePipeline2/Programs</u> folder.
- [bwa](http://bio-bwa.sourceforge.net/)
- [picard](https://broadinstitute.github.io/picard/)
- [fastp](https://github.com/OpenGene/fastp)
- [seqtk](https://github.com/lh3/seqtk)
- [kraken](http://ccb.jhu.edu/software/kraken/)
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- [snp-sites](https://github.com/sanger-pathogens/snp-sites)
- [pgzip](https://zlib.net/pigz/)
- [MultiQC](https://multiqc.info/)
- [QualiMap](http://qualimap.conesalab.org/)
- [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/4418062771227-Mutect2)
- [SnpEff](http://pcingola.github.io/SnpEff/)

Python libraries:
- [pandas](https://pandas.pydata.org/)
- [PyVCF](https://pyvcf.readthedocs.io/en/latest/)
- [pairsnp](https://github.com/gtonkinhill/pairsnp)

### Structure

ThePipeline2 folder is structured as follows:

- ThePipeline_v2
    - data
      - Configs
      - libs
      - Paths
    - Programs
      - bedtools2
      - bwa
      - fastp
      - gatk-4.2.5.0
      - Kraken
      - MultiQC_env
      - picard
      - pigz
      - Pilon (not in use)
      - qualimap_V2.2.1
      - samtools-1.15
      - seqtk
      - snpEff
      - snp-sites
    - PipeModules

#### data

This folder contains data needed for the different modules to work properly. 

You can find the **catalogWHO2021_antibiotics_customized.tsv** file (needed for the [resistance](#detect-genomic-determinants-of-resistance-resistance) module), **H37Rv.annotation_new.tsv** file (needed for the [annotation_filter](#filter-by-annotation-annotation_filter) module) and **snp_phylo.tsv** and **snp_phylo_reduced.tsv** files (needed for the [typing](#lineage-typing-typing) module). 

In the *Configs* folder, you will find the [fastp](https://github.com/OpenGene/fastp) tool configuration file (**fastp_default.config**), the [MultiQC](https://multiqc.info/) tool configuration file (**multiqc.config**), and information about the software versions used (**software_versions.txt**). The software_versions.txt file is used to write the PREFIX.history files, and useful to track program versions used to analyze the different datasets. 

In the *libs* folder there is a version of the libstdc++ library (**libstdc++.so.6**), to avoid problems with the [coverage](#coverage-calculation-and-filtering-coverage) module that sometimes arise for the presence of old versions of this library in the servers. You can find also the **pairsnp-python**, cloned from [the pairsnp GitHug repository](https://github.com/gtonkinhill/pairsnp), and necessary for the [distances](#calculate-genetic-distances-distances) module to run.


In the *Paths* folder, you will find two different files; **programs_path** in which the relative path to each of the programs used by ThePipeline2 is specified, and **data_path** which contains the absolute path to the MTBC ancestor reference, the kraken database, and the picard dictionary.  


#### Programs
In this folder you you will find all the programs needed by ThePipeline2 to work. The folders containing the programs BWA, seqtk and snp-sites are zipped in the GitLab repository, so you will need to unzip them after downloading the repository and before running ThePipeline2.

#### PipeModules
This folder contains the python code of each of ThePipeline2 modules.


### Instalation

ThePipeline2 has been written in python 3. So, you will need to have python 3 installed in your system to run ThePipeline2. This pipeline is prepared to work in the Tuberculosis Genomic Unit HPC servers. If you intend to run this pipeline in other machines, you will need to:
1. Download ThePipeline2 repository from [GitLab](https://gitlab.com/tbgenomicsunit/ThePipeline)
2. Check that precompiled binaries of the programs (/Programs folder) work in your specific computer. The folders containing the programs BWA, seqtk and snp-sites are zipped in the GitLab repository, so you will need to unzip them at this stage.
3. Change the absolute paths of the /data/Paths/data_path file to point to the kraken database, the MTBC ancestor genome and the picard dictionary. The MTBC most recent ancestor inferred genome can be found [here](http://tgu.ibv.csic.es/wp-content/uploads/2019/09/MTB_ancestor_reference.fasta.gz), and cited from [Comas et al., 2010](https://www.ncbi.nlm.nih.gov/pubmed/20495566/).
4. You may need to install the [pairsnp python library](https://github.com/gtonkinhill/pairsnp) for the [ThePipeline2 distances](#calculate-genetic-distances-distances) module to work. Follow the instructions in /ThePipeline_v2/data/libs/pairsnp-python/README.md
5. You may need to reinstall the MultiQC envinronment for [ThePipeline2 multiqc](#generate-multiqc-report-multiqc) module to work. This could be easily achived by executing the following in /ThePipeline_v2/Programs/ folder:
    ```bash
    rm -rf MultiQC_env/
    python3.7 -m venv MultiQC_env
    source MultiQC_env/bin/activate 
    pip install multiqc

    # sometimes installation fails showing an error with the Pillow package
    # in this case, execute
    pip install --upgrade pip
    pip install --upgrade Pillow

    # try again to install multiqc
    pip install multiqc
    ```
6. You may need to install [pandas](https://pandas.pydata.org/) and [PyVCF](https://pyvcf.readthedocs.io/en/latest/) libraries for [ThePipeline2 calling](#variant-calling-calling) module to work. This could be easily achived by executing the following:
     ```bash
    pip3 install --user pyvcf
    pip3 install --user pandas
    ```

Sudo privilegies are desirable to install ThePipeline2, although it may be installed in local folder without any problem (not tested).

### Updating ThePipeline2

It may happen that, at some point, some programs have to be updated. Please, bear in mind that this step is critic and may compromise ThePipeline2 functionality. If after careful evaluation of pros- and cons- you decided to go ahead with the update, follow these steps:
1. First, check that the output of the new version of the program has exactly the same format of the previous one. ThePipeline2 has been written to work with specific types of data and formats. Changes in the formats could lead to execution errors. Or worst, it may lead to unnotified errors or malfunctions not detectable at first sight.
2. Create a new folder in /Programs folder with the new version of the program. I recommend to leave the old version in a 'legacy_programs' folder, just in case you need to return to it. 
3. Change the **/data/Paths/programs_path** file, so the program acronym now points to the new version of the executable. For example, if you update the gatk version to 4.3, you must change the gatk line to something like *gatk=gatk=/Programs/gatk-4.3.0/gatk*.
4. Change the **/Config/software_versions.txt** file to match the new version of the program. This info is used to write the PREFIX.history files, and useful to track program versions used to analyze the different datasets. 
5. If you update [fastp](https://github.com/OpenGene/fastp) or [MultiQC](https://multiqc.info/), you would need to update their config files (**/Config/fastp_default.config** and **/Config/multiqc.config**).
6. If the arguments of the program do not match those used by the previous version, you will need to go to the /PipeModules folder and change the python code of the modules affected. Bellow, in each module section, you have a description of the programs used on them. This is a critical step and should be executed by a person keen on python programming. The calls to external programs are executed using the *subprocess.call* function. Arguments are passed in the form of a list of strings, usually called **cmd_xxx**, in which each element of the list is an argument. You will need to modify this list. Do not use a single string with arguments separated by spaces! It will make the method function to crash. Following the previous example, imagine that gatk version 4.3 has changed the way in which Mutect2 define the number of threads to use, and the argument previously called "--native-pair-hmm-threads" now is called '--threads'. You must go to CallingMutect2.py and change the following lines:
   ```python
   # old version
   cmd_Mutect2 = [gatk, "Mutect2",
                   "-R", reference, "-I", "{}{}".format(prefix, ext),
                   "-O", "{}.vcf".format(prefix), "-OVI", "false",
                   "--verbosity", "ERROR",
                   "--QUIET", "true", "-mbq", min_qual,
                   "--callable-depth", min_depth,
                   "-mnp-dist", "0",  # important to not filter MNPs later
                   "--linked-de-bruijn-graph", "true",
                   "--native-pair-hmm-threads", threads]
    # new version
    cmd_Mutect2 = [gatk, "Mutect2",
                   "-R", reference, "-I", "{}{}".format(prefix, ext),
                   "-O", "{}.vcf".format(prefix), "-OVI", "false",
                   "--verbosity", "ERROR",
                   "--QUIET", "true", "-mbq", min_qual,
                   "--callable-depth", min_depth,
                   "-mnp-dist", "0",  # important to not filter MNPs later
                   "--linked-de-bruijn-graph", "true",
                   "--threads", threads]
   ```
   If you neew to add new arguments instead of change current ones, just add as elements to the list as new arguments needed.
 
---


# Running ThePipeline2

At this moment (07/2022) we have three different High Performance Servers at the TGU: Koch, Yersin and Salas. All of them have enough RAM memory and CPU capacity for running all _ThePipeline2_ modules. However, take into account that [kraken](http://ccb.jhu.edu/software/kraken/) needs hundreds of GB of RAM for running, so if you intend to run ThePipeline2 in other machine the [ThePipeline2 kraken](#taxonomic-filtering-kraken) module may not work if you not meet a minimum of [174GB RAM](https://ccb.jhu.edu/software/kraken/MANUAL.html#system-requirements). 

'Standard' runs consist on the execution of the following sequential steps:

1. [FASTQ Filtering (ThePipeline2 fastclean)](#fastq-filtering-fastclean)
2. [Taxonomic filtering (ThePipeline2 kraken)](#taxonomic-filtering-kraken)
3. [Mapping (ThePipeline2 mapping)](#mapping-mapping)
4. [Variant calling (ThePipeline2 calling)](#variant-calling-calling)

Each of these steps uses as an input the output of the previous one. Other modules complement the 'core' ones:

- [Coverage calculation and filtering (ThePipeline2 coverage)](#coverage-calculation-and-filtering-coverage) (depends on the output generated by ThePipeline2 mapping)
- [Filter by annotation (ThePipeline2 annotation_filter)](#filter-by-annotation-annotation_filter) (depends on the output generated by ThePipeline2 calling)
- [Generate multifastas (ThePipeline2 consensus)](#generate-multifastas-consensus) (depends on the output generated by ThePipeline2 calling and ThePipeline2 annotation_filter)
- [Organizing results (ThePipeline2 organize)](#organizing-results-organize)
- [Detect genomic determinants of resistance (ThePipeline2 resistance)](#detect-genomic-determinants-of-resistance-resistance) (depends on the output generated by ThePipeline2 calling)
- [Generate MultiQC report (ThePipeline2 multiqc)](#generate-multiqc-report-multiqc)(depends on the output generated by ThePipeline2 fastclean, ThePipeline2 kraken and ThePipeline2 mapping)
- [Lineage typing (ThePipeline2 typing)](#lineage-typing-typing) (depends on the output generated by ThePipeline2 calling)
- [Calculate genetic distances (ThePipeline2 distances)](#calculate-genetic-distances-distances) (depends on the output generated by ThePipeline2 consensus)
- [Calculate transmission clusters (ThePipeline2 getclusters)](#calculate-transmission-clusters-getclusters) (depends on the output generated by ThePipeline2 distances)

Some of the modules must be executed individually over the samples (fastclean, kraken, mapping, calling, coverage annotation_filter) while some others must be executed over a complete set of samples indicating the folder in which the data is placed (consensus, organize, resistance, multiqc). See the manual section of each module for specific details on how to run each module.

For those modules that must be executed individually on each sample, the best option would be the use of the [xargs](https://ss64.com/bash/xargs.html) bash tool to paralelize the run. It will work with something like:

```bash
ls *.cram | cut -f1 -d'.' | xargs -I sample -P 4 sh -c "ThePipeline2 calling -p sample -t 5"
```

In this case, you listed all the CRAM files, keep the first part of the file name (before the first '.') and xargs save it in the 'sample' variable. Later, it executes the code between quotes. The -P option can be used to parallelize the operation, so the code between brackets will be executes P times (4 in this case), and each time the 'sample' variable will has a different value from the complete list file names. 

Bear in mind that the maximum numer of processes/threads dedicated to each run must be in accordance with your system characteristics. In this particular case, you will use 20 CPU nodes (5 threads x 4 processes). The use of the [nice](https://ss64.com/bash/nice.html) instruction is highly recommended!!
```bash
ls *.cram | cut -f1 -d'.' | xargs -I sample -P 4 sh -c "nice -n 5 ThePipeline2 calling -p sample -t 5"
```

Finally, ThePipeline2 generates, for each input data, a .history file in which all ThePipeline2 processes executed, including date and time, programs, programs version and arguments are detailed. 

---

# Modules
## FASTQ Filtering (fastclean)
This module relies on the [fastp](https://github.com/OpenGene/fastp) tool for an initial preprocessing of the fastq files. It is run by typing

```bash
ThePipeline2 fastclean
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -f, --fastq| [FASTQ [FASTQ ...]]|The input fastq files. They can be compressed in the fastq.gz form|
| -c, --config| CFGFILE |plain text file with fastp parameters to be used. By default: --cut_tail, --cut_window_size=10, --cut_mean_quality=20, --length_required=50, --correction. Default config file can be found in the config folder. |
| -t, --threads | THREADS | Number of threads to be used. Default: 1|
| -v, --verbose ||Print the command line used to call fastp in terminal. Deafult: false|
| --phread64 ||Use Phred64 scale for input fastq. Output will be in Phred33|
| -p, --prefix |PREFIX| Prefix/name of the sample. It will be used for naming the outputs generated|


This module takes as an input a fastq file(s) and returns another fastq file(s) with the reads that passed the filters defined in the config file. By default, the following filters are applied:

- --cut_tail : move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise
- --cut_window_size = **10** : the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000.
- --cut_mean_quality = **20**  : the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36
- --length_required = **50** : reads shorter than length_required will be discarded
- --correction : enable base correction in overlapped regions (only for PE data). If an proper overlap is found, it can correct mismatched base pairs in overlapped regions of paired end reads, if one base is with high quality while the other is with ultra low quality. If a base is corrected, the quality of its paired base will be assigned to it so that they will share the same quality.

There are many other filters that can be applied, just check the [fastp](https://github.com/OpenGene/fastp) manual, create a new config file and pass it to ThePipeline fastclean with -c/--config option. If everything was fine, the *fastclean* module generates the following output:
- PREFIX.(P1/2).clean.fastq.gz : FASTQ file containing the reads that passed all filters. In the case of single-end data, only 1 file is generated. When dealing with paired-end data, each of the fastq files will be labeled as P1 or P2.
- PREFIX.cleanlog: A text file summarizing the filters applied and the reads that passed each of the filters. Will be stored in the /reports folder
- PREFIX.fastp.json: A json file, with information about the filters applied and the reads that passed each of the filters. Will be stored in the /reports folder 

**Examples of usage:**

Basic
```bash
ThePipeline2 fastclean -f ERR9029927_[12].fastq.gz -p ERR9029927 -t 4 -v
```

Redefining filtering options, saving to a config file and passing it to the program
```bash
echo -e "--cut_tail\n--cut_window_size=10\n--cut_mean_quality=30\n--length_required=100\n--correction \n" > fastp_config.txt;

ThePipeline2 fastclean -f ERR9029927_[12].fastq.gz -p ERR9029927 -t 4 -v -c fastp_config.txt
```

---


## Taxonomic filtering (kraken)
This module relies on [kraken](http://ccb.jhu.edu/software/kraken/), [seqtk](https://github.com/lh3/seqtk) and [pgzip](https://zlib.net/pigz/). It is used for taxonomic classification and/or filtering of the reads contained in a fastq file. 

```bash
ThePipeline2 kraken
```

It has three main modes of running, that can be runned independently or at the same time. If you want to run them independetly, bear in mind that **filter** and **report** use as an input the results derived from **classify**:
- classify : Is the first step. Kraken assigns a taxonomic level to each read in the fastq files by comparing to a previously constructed database. Results are written to .kraken files. 
- filter : Uses the .kraken files generated by the **classify** module for filtering reads according to the taxonomic name given as an argument, using grep and seqtk. These reads are stored in a new 'filtered' fastq file, compressed using pgzip.  
- report : Uses the .kraken files generated by the **classify** module to generate several human-readable report files, in tabular format.

Main arguments:

|Argument|Value|Description|
|:------------- |:---|:-----|
| --classify ||Run kraken to classify fastq reads using krakenDB|
| -f, --fastq| [FASTQ [FASTQ ...]]|The input fastq files. They can be compressed in the fastq.gz form. It is a compulsory argument for the '--classify' and '--filter' options.|
| -c, --compressed|  | Whether the input fastq files are gzip compressed.  |
| --paired|  | Whether the input fastq files are paired. |
| --db| PATH | Kraken database to use. By default uses the one stored in /data/Databases/KrakenDB/KRefSeqHuman |
| --preload|  | Load the kraken database in RAM memory. |
| -t, --threads | THREADS | Number of threads to be used. Default: 1|
| --filter ||Save reads that matched the -m/--matching argument into .filtered.fastq.gz|
| -m, --matching| NAME | Taxonomic name for matching. By default use _Mycobacterium tuberculosis_.|
| -r, --report || Generates two report files summarizing the results obtained by kraken. One file summarizes the potential contaminants at the genus level while other summarizes the contaminants at the species level|
| -p, --prefix |PREFIX| Prefix/name of the sample. It will be used for naming the outputs generated|
| --noclean ||Do not remove .labels, .readlist and .kraken files. By default, these files are removed at the end of each execution.|

If everything was fine, the kraken module generates the following output:
- PREFIX.kraken (removed by default, obtained with --classify): Default kraken output. Contain five tab-delimited fields; from left to right, they are:
    -   "C"/"U": one letter code indicating that the sequence was either classified or unclassified.
    -   The sequence ID, obtained from the FASTA/FASTQ header.
    -   The taxonomy ID Kraken used to label the sequence; this is 0 if the sequence is unclassified.
    -   The length of the sequence in bp.
    -   A space-delimited list indicating the LCA mapping of each k-mer in the sequence. For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
        -   the first 13 k-mers mapped to taxonomy ID #562
        -   the next 4 k-mers mapped to taxonomy ID #561
        -   the next 31 k-mers contained an ambiguous nucleotide
        -   the next k-mer was not in the database
        -   the last 3 k-mers mapped to taxonomy ID #562
- PREFIX.labels (removed by default, obtained with --classify): Default kraken-translate output. Is a text file with two tab-delimited columns, and one line for each classified sequence; unclassified sequences are not reported. The first column are the sequence IDs of the classified sequences, and the second column contains the taxonomy of the sequence.
- PREFIX.filtered.readlist (removed by default, obtained with --classify): Contains the ID of the reads in PREFIX.labels whose taxonomical classiffication matches the -m/--matching argument (by default _Mycobacterium tuberculosis_)
- PREFIX.kraken.report (obtained with --classify): Default kraken-report output. This file is tab-delimited, with one line per taxon. The fields of the file, from left-to-right, are as follows:
    - Percentage of reads covered by the clade rooted at this taxon
    - Number of reads covered by the clade rooted at this taxon
    - Number of reads assigned directly to this taxon
    - A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
    - NCBI taxonomy ID
    - indented scientific name 
- PREFIX.P[12].filtered.fastq.gz (obtained with --filter): Compressed FASTQ files containing the reads whose taxonomical classification matches the -m/--matching argument (by default _Mycobacterium tuberculosis_)
- PREFIX.nfilter (obtained with --filter): Text file containing a short report with the total number of reads analyzed on the original FASTQ file and how many of them belong to the MTBC.
- PREFIX.genus.contaminants (obtained with --report): Tab-delimited file, formatted from the PREFIX.kraken.report file. Contains information about the different taxonomic elements detected in the original FASTQ file, at the genus level.
- PREFIX.species.contaminants (obtained with --report): Tab-delimited file, formatted from the PREFIX.kraken.report file. Contains information about the different taxonomic elements detected in the original FASTQ file, at the species level.

**Examples of usage:**

Classify
```bash
ThePipeline2 kraken --classify -f ERR9029927.P[12].clean.fastq.gz --paired --compressed -p ERR9029927 -t 20 --noclean
```

Filter
```bash
ThePipeline2 kraken --filter -p ERR9029927 --noclean
```

Report
```bash
ThePipeline2 kraken --report -p ERR9029927
```

Complete
```bash
ThePipeline2 kraken -f ERR9029927.P[12].clean.fastq.gz --paired --compressed -p ERR9029927 --report --classify --filter -t 20
```

Changing the taxonomic name for matching
```bash
ThePipeline2 kraken -f ERR9029927.P[12].clean.fastq.gz --paired --compressed -p ERR9029927 --report --classify --filter -t 20 -m "Mycobacterium smegmatis"
```

To speed up the execution of the --classify module it is advisory to load the kraken database in the RAM memory. This can be achieved with the --preload argument. Use the preload argument only in the first analysis, and the database will remain in RAM for the following executions.
```bash
ThePipeline2 kraken --classify --filter --report -f ERR9029927.P[12].clean.fastq.gz --paired --compressed -p ERR9029927 -t 20 --preload;

ThePipeline2 kraken --classify --filter --report -f ERR9029928.P[12].clean.fastq.gz --paired --compressed -p ERR9029927 -t 20
```
---

## Mapping (mapping)
This module relies on [bwa](http://bio-bwa.sourceforge.net/), [picard](https://broadinstitute.github.io/picard/), [samtools](http://www.htslib.org/) and [QualiMap](http://qualimap.conesalab.org/). Input fastqs are mapped to a reference genome with bwa. The resulting aligment is sorted, and optical and PCR duplicates are removed with picard. A quality report is generated with QualiMap. Finally, the aligment is stored into a [CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format)#:~:text=CRAM%20is%20a%20compressed%20columnar,Map%20(BAM)%20file%20formats) file (compressed BAM).
```bash
ThePipeline2 mapping
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -f, --fastq| [FASTQ [FASTQ ...]]|The input fastq files. They can be compressed in the fastq.gz form|
| -r, --reference| REFERENCE | Reference genome in FASTA format. By default, the MTBC inferred ancestor. If other reference genome is selected, take into account that you will first need to index the genome using BWA and GATK|
| -p, --prefix |PREFIX| Prefix/name of the sample. It will be used for naming the outputs generated|
| -t, --threads | THREADS | Number of threads to be used. Default: 1|
| -i, --index ||Index reference genome (only first time) using the 'bwtsw' algorithm.|
| --keep-dupbam||Keep the BAM file with the duplicates not removed.|
| --no-dedup ||Do not remove duplicates. Duplicates are removed by default using picard-tools.|
|-mapq|MAPQ_CUTOFF|Mapping quality filter. Reads with MAPQ below this value will be filtered out in the BAM file. DEFAULT=60|

If everything was fine, the mapping module generates the following output:
- PREFIX.sort.cram: Sorted CRAM file (compressed BAM) with the fastq reads mapped to the reference genome. By default, this file has the duplicated reads detected removed. If mapping has been runned with the --no-dedup option, the potential duplicated reads will be present.  
- PREFIX.dup.metrix: Report of the read duplication removal process. It shows information about the number of optical and PCR duplicates detected and removed by picard-tools.
- PREFIX.dup.sort.bam (obtained with --keep-dupbam): 'Original' BAM file, after the duplicated reads have been removed.
- PREFIX.sort.bam folder: This folder is generated with the QualiMap tool. It contains information about the BAM quality parameters. This folder will be later used by the multiqc module.

Please, bear in mind that if the --no-dedup argument is present, you will obtain only one output aligment (PREFIX.sort.cram) with the duplicated reads on it. On the contrary, if you select --keep-dupbam you will obtain two output alignments, one with the duplicated reads removed (PREFIX.sort.cram) and other with the duplicated reads present (PREFIX.dup.sort.bam).
Also, take into account that QualiMap needs a graphic envinronment to work properly. Hence, you will need to connect to the servers using the **-X argument in the ssh instruction** when running ThePipeline2.

**Examples of usage:**

Basic
```bash
ThePipeline2 mapping -f ERR9029927.P[12].filtered.fastq.gz -t 4 -p ERR9029927
```

Requiring a lower MAPQ value of 40
```bash
ThePipeline2 mapping -f ERR9029927.P[12].filtered.fastq.gz -t 4 -p ERR9029927 -mapq 40
```

When using a new reference, you would need to index the FASTA file before first usage
```bash
ThePipeline2 mapping -i -r ~/H37Rv.fasta -f ERR9029927.P[12].filtered.fastq.gz -t 20 -p ERR9029927
```
---

## Variant calling (calling)

This module mainly relies on [gatk](https://gatk.broadinstitute.org/hc/en-us/articles/4418051394587--Tool-Documentation-Index). It uses the Mutect2, FilterMutectCalls and SelectVariants tools to perform the variant calling step. It works as follows:
- Mutect2 generates a single VCF file with all variation found.
- FilterMutectCalls adds filtering information to the VCF generated in the previous step.
- This VCF is later splitted in two files: one with all the INDELS detected (saved in INDEL.VCF file), and one with the variants identified by GATK as not being INDELS, MNPs (phased SNPs), SYMBOLIC or NO_VARIATION (saved in the SNP.VCF file). The original VCF file is deleted by default, althoug it can be stored for debugging purpouses with the corresponding argument (-k).
- The remaining VCF files are later filtered and reformated:
  - all the mutiallelic positions found in the SNP.VCF file (if any) are moved to a new file called MULTIALLELIC.VCF. SNPs found in multiallelic positions remain also in the SNP.VCF file, but only with the information about the SNP with higher frequency. 
  - the SNP.VCF file is modified and all the positions having AF <= 'filtering frequency cut-off' will have an 'hetero' code in the ALT field, in an attempt to mimic the VarScan behaviour in the previous version of ThePipeline.
  - SNPs above frequency and depth cut-offs are stored in a new VCF called EPI.SNP.VCF file.
  - SNPs under a minimum frequency and a minimum depth are removed from the SNP.VCF file.
  - INDELS under 10 depth and under minimum frequency are removed from the INDEL.VCF file.
- Later, we applied an extra filtering step, in which we removed SNPs in EPI.SNP.VCF file that accumulate in small regions. Using a sliding window, we evaluate the density of SNPs and remove those above a threshold.
- Later, using bedtools we infer the low coverage positions (those having less than 3X depth).
- Next, we inferred the WT positions as those that are not lowcov, not present in the SNP.VCF, not present in the INDEL.VCF file and not present in the MULTIALLELIC.VCF file. 
- Finally, the VCF files (except MULTIALLELIC.VCF) are annotated using [SnpEff](http://pcingola.github.io/SnpEff/). By default, SnpEff is executed automatically when the reference genome used is the MTBC inferred ancestor. If other reference is passed by the -r argument, SnpEff is not automatically executed, hence users will need to annotate VCF files 'manually' (you can use SnpEff too, downloading or generating a new database for the reference genome used).

![calling-schema](calling.png)

Please, read bellow for a more detailed description of the files derived.

ATTENTION: If you tried to run the module in a folder in which the calling module has been previously run (meaning that there are .snp.vcf, indel.vcf and calling-derived files) it may crash due to conflict with the output file names in some of the steps. If you want to re-run the calling module, first move or delete the outputs of the previous runs.  

```bash
ThePipeline2 calling
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -r, --reference| REFERENCE | Path to the reference genome in FASTA format. By default, the MTBC inferred ancestor.|
| -p, --prefix | PREFIX | Prefix/name of the sample. It will be used for naming the outputs generated|
| -e, --extension | {.sort.bam,.sort.cram} | Whether the input files are in bam or cram format. By default: .sort.cram|
| -t, --threads | THREADS | Number of threads to be used. Default: 1|
| -gVCF, --genomic_VCF |  | Generates a genome-wide VCF and then the calling process stops. For debuging purpouses only. |
| -k, --keep_vcf | | Keeps the original VCF generated by Mutect2. By default this file is removed. |
| -min_d, --min_depth | MIN_DEPTH | Minimum depth of coverage to consider a position callable. By default 3 reads. |
| -min_q, --min_qual | MIN_QUAL | Minimum base quality to consider a position callable. By default 15. |
| -min_f, --min_freq | MIN_FREQ | Minimum frequency to consider an alternative allele. By default 0.05. |
| -filt_d, --filter_depth | FILTER_DEPTH | Minimum depth to include a SNP in the EPI.snp file. By default 20 reads. |
| -filt_f, --filter_freq | FILTER_FREQ | Minimum alternative allele frequency to include a SNP in the EPI.snp file. By default 0.9 |
|--skip_dens_filt | | Skip the density filter. By default is activated |
|-w, --window | WINDOW | Size of the slidding window used for the density of SNPs filtering. By default 10 pb.|
|-dens, --density | DENSITY | Density cut-off. If more or equal number of SNPs than the cut-off are found in the defined slidding window, they will be removed. By default 3.|

If everything was fine, the calling module generates the following output:
- PREFIX.gvcf (not generated by default): This file will be only generated by the use of the -gVCF/--genomic_VCF option. It is a genome-wide VCF, with all the variation detected in every reference position, before applying any filtering.
- PREFIX.vcf (not keep by default): This file is the original VCF file generated by Mutect2. By default it is removed, but with the option -k/--keep_vcf the file remains. 
- PREFIX.lowcov.tsv: Text file containing the bases that didn't meet the minimum read depth for being considered as callable (defined by -min_d, by default 3X).
- PREFIX.snp.vcf: VCF file containing all bases having an alternative nucleotide in comparison with the reference genome. This file is annotated. In this file you will find the SNPs detected by Mutect2 that:
  - Are above min_f and min_d
  - SNPs with freqs between min_f and filt_f, will have an undetermined position in the ALT field, similar to what happened in the past with ThePipeline (v.1)
  - REF and ALT columns must be 1 nucleotide long. In multiallelic positions having insertions, only the SNP allele (if present) is included in this file.
- PREFIX.indel.vcf: VCF file containing small indels. This file is annotated. Indels in this file are those detected by mutect having depth above 10X and freq above min_f.
- PREFIX.multiallelic.vcf: All positions found by Mutect2 that are multiallelic. Some positions may appear in both: PREFIX.multiallelic.vcf and PREFIX.snp.vcf. BUT, in the PREFIX.snp.vcf will only appear information about the SNPs, and a message stating MULTIALLELIC in the INFO field of the VCF file. 
- PREFIX.wt.tsv: Text file containing all the positions that are above mind_d, and that are not present in PREFIX.lowcov.tsv, PREFIX.snp.vcf, PREFIX.indel.vcf and PREFIX.multiallelic.vcf.
- PREFIX.EPI.snp.vcf: VCF annotated file containing all the SNPs from the PREFIX.snp.vcf file that passed the following filters:
  - PASS in the FILTER field of the PREFIX.snp.vcf file (flag added by FilterMutectCalls) 
  - Depth equal or above -filt_d (20 by default).
  - Frequency equal or above -filt_f (0.9 by default).
  - Not in high density regions (unless density filter had been deactivated).
- PREFIX.dens_removed_snps.vcf: VCF file containing all the SNPs from the PREFIX.EPI.snp.vcf file (meaning that they have passed depth and freq thresholds) that have been removed by the density filter. By default, remove SNP if -dens (3 by default) or more SNPs are found in -window (10 by default) pair bases. This file is annotated.   

**Examples of usage:**

Basic

```bash
ThePipeline2 calling -p ERR9029927 -t 20
```

Relax filtering cut-offs. Filter out from the PREFIX.EPI.snp.vcf file those SNPs having less than 10X and under 0.1 freq. Do not apply density test 

```bash
ThePipeline2 calling -p ERR9029927 -t 20 -filt_d 10 -filt_f 0.1 --skip_dens_filt
```


---
## Coverage calculation and filtering (coverage)
This module uses [bedtools](https://bedtools.readthedocs.io/en/latest/) to calculate the sequencing depth at each genomic position of the reference genome. Later, using this information the module calculates coverage and depth. It can also filter out samples that do not reach minimum coverage and median/mean depths cut-offs. 
```bash
ThePipeline2 coverage
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -r, --reference| REFERENCE | Path to the reference genome in FASTA format. By default, the MTBC inferred ancestor. |
| -p, --prefix |PREFIX| Prefix/name of the sample. It will be used for naming the outputs generated|
| -e, --extension | {bam,cram} | Whether the input files are in bam or cram format. Default: cram|
| -f, --filter ||Activate the coverage filter. Samples not passing thresholds are moved to NoPassCov/|
| --min-depth-mean |VALUE|Minimum mean depth for considering a sample. Default: 0|
| --min-depth-median |VALUE|Minimum median depth for considering a sample. Default: 20|
| --min-coverage |VALUE|Minimum genomic coverage for considering a sample. It must be a value in range [0.0, 1.0]. Default: 0.95|
| --keep-coverage||Keep .coverage files|
| --deepth4cov | VALUE |Minimum depth to consider covered a base. Default: 10|

If everything was fine, the coverage module generates the following output:
- PREFIX.coverage (removed by default): It is a tab-delimited file obtained after running the genomeCoverageBed module of the bedtools2 package. It has 3 columns that correspond to:
    1. Chromosome
    2. Genomic position in the chromosome
    3. Depth of features (reads) overlapping at this chromosome position
    
    By default this file is removed, but its removal can be skipped by using the --keep-coverage option.
- PREFIX.meancov: A tab-delimited file with the three statistics calculated by the module: mean depth, median depth and genomic coverage.
- NoPassCov/ (obtained with -f/--filter): All the files derived from samples that do not pass the --min-depth-mean, --min-depth-median or --min-coverage cut-offs will be moved to this folder. If all samples pass filters, this folder will not be created.

**Examples of usage:**

Basic
```bash
ThePipeline2 coverage -p ERR9029927 -e cram
```

Filtering samples with genomic coverage under 90%
```bash
ThePipeline2 coverage -p ERR9029927 -e cram -f --min-coverage 0.90
```
---

## Filter by annotation (annotation_filter)
This module filter VCF files, removing positions called in problematic regions (mainly PE/PPE/repeats/phages). These problematic regions are marked as *DISCARD* in the file **H37Rv.annotation_new.tsv** that can be found in the ThePipeline2/data folder. You can download this file and find more information about how these regions to mask were derived in the [TGU repository webpage](http://tgu.ibv.csic.es/?page_id=1794), in the section 'Genomic regions masked in phylogenetic and epidemiological studies'. In principle, all the VCF files generated by the [calling](#variant-calling-calling) module can be filtered, as the **annotation_filter** module just evaluates the column 'POS' of the VCF files. 
```bash
ThePipeline2 annotation_filter
```

This module only has one argument
|Argument|Value|Description|
|:------------- |:---|:-----|
| -s, --snp-file| VCF | VCF file to be filtered |

If everything was fine, the annotation_filter module generates the following output:

- VCF.annoF: VCF file with the positions falling in the problematic regions indicated in the **H37Rv.annotation_new.tsv** file removed.

**Examples of usage:**

Basic
```bash
ThePipeline2 annotation_filter -s ERR9029927.EPI.snp.vcf
```

---

## Generate multifastas (consensus)
This module generates, for each sample, a fasta file with the variants called in the [calling](#variant-calling-calling) module. It accepts multiple folders as an input. The module scans all these folders and uses the EPI.snp.vcf files found on them to construct the multifastas, following these steps:
  
  1. Generates a SNP_table, joining all the SNPs contained in all the EPI.snp.vcf files found in the input folders. The SNP_table contains a non redundant annotated list of all these polymorphisms. Later, for each sample, generates a FASTA containing all the SNPs called in this sample by doing:
     - Initiallize the FASTA file, creating a sequence containing the REF alleles of the previously generated SNP_table. Later, for each position in sequence:  
        1. If the position is in the **wt.txt** file, leave the reference allele, otherwise continue.
        2. If the position is in **.snp.vcf** file, put the ALT allele, otherwise continue.
        3. If the position is in **.indel.vcf**, put '-', otherwise continue.
        4. If none of the above, put an 'N'
  2. Joins all individual FASTA a single multifasta file. Remove the individual FASTA files.
  4. Uses snp-sites to remove the multifasta conserved positions and generates a new 'snp-sites multifasta' with its associated 'snp-sites SNP_table'
  5. Removes the resistant-confering positions (present in the WHO catalogue), generating a new 'snp-sites no-resis' multifasta and its associated 'snp-sites no-resis SNP_table'

```bash
ThePipeline2 consensus
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -i | PATHS | Paths to the folders containing all the samples that have to be included in the complete multifasta. Default: current folder |
| -t, --threads | THREADS | Number of threads to use. Default = 1.|
| -p, --prefix |PREFIX| It will be used for naming the outputs generated. By default current date (YearMonthDay)|

If everything was fine, the annotation_filter module generates the following output:

- PREFIX.mf.fasta: Initial multifasta file
- PREFIX.mf_gap.fasta: From the initial multifasta file, all the non-ATGC bases are substituted by gap symbols '-'.
- PREFIX.SNP_table.txt: Non redundant list of all the SNPs that passed all filters in, at least, one sample of the analyzed ones (mean that they are present in at least one EPI.snp.vcf). Each row correspond to a position in the PREFIX.mf.fasta file (except for the header). It also includes information about the gene in which the SNP has been called, and annotation information. In the case of 'multiallelic' positions in the multifasta, all the alternative alleles are reported in the ALT field, as well as the alternative annotations. If one genomic position belongs to more than one genomic region (overlapping genes for example), in the 'Rv_number' will appear the flag '+OL'. The annotation information of the SNP, in this case, will reffer to the first feauture of the overlapping regions.
- PREFIX.mf_gap.snp-sites.fasta: FASTA file derived from PREFIX.mf_gap.fasta file, only containig variable positions.
- PREFIX.SNP_table.snp-sites.txt: SNP_table corresponding to the PREFIX.mf.snp-sites.fasta file.
- PREFIX.mf_gap.snp-sites.no-resis.fasta: FASTA file derived from PREFIX.mf_gap.snp-sites.fasta file, only containig variable positions and without the SNP resistant positions that appear in the [WHO catalogue of mutations](https://www.who.int/publications/i/item/9789240028173).
- PREFIX.SNP_table.snp-sites.no-resis.txt: SNP_table corresponding to the PREFIX.mf.snp-sites.no-resis.fasta file.
- PREFIX_invariants.txt: Text-file containig the MTBC genome invariant sites (A, C, G, T) corresponding to the PREFIX.mf.snp-sites.no-resis.fasta alignment. This file is intended to be used with IQTREE and the option -fconst. It can be passed to IQTREE directly: 

```bash
iqtree -s PREFIX.mf_gap.snp-sites.no-resis.fasta -m GTR **-fconst $(grep -v '#' PREFIX_invariants.txt)** -bb 1000
```

**Examples of usage:**

Basic
```bash
ThePipeline2 consensus -t 15
```

Including other folder
```bash
ThePipeline2 consensus -t 15 -i ./ /data2/reference_dataset/representative_genomes_MTBC_Sebastien/zSNPs/
```

---

## Organizing results (organize)
This module organizes most of the outputs obtained by the different _ThePipeline2_ modules in several folders.
```bash
ThePipeline2 organize
```

---

## Detect genomic determinants of resistance (resistance)
This module compares the different variants identified by the [calling](#variant-calling-calling) module againts a catalog of known drug-resistance (DR) conferring mutations ([WHO catalog](https://www.who.int/publications/i/item/9789240028173)). It scans all the samples present in the folder in which it is executed. It returns a report with the potential resistant profile of each sample and the genetic determinants identified. It works as follows. For each sample:
1. List all the potential genetic derminants (only SNPs and Dels at this moment, it will also evaluate double SNPs on the same codon) in the catalog.
2. If all the catalog positions are in the .wt.txt file, write the potential resistant profile. If not, 3. 
3. Check if the remaining positions are in the lowcov.vcf file, and write the potential resistant profile. If not, 4
4. Check if the remaining positions are in the .snp.vcf file or in the indel.vcf file. Write the resistant profile of the sample and a report with information about the DR conferring positions found.

It has an alternative mode of running. If the '-ns' option is selected, it looks only for SNPs falling in the DR genes listed in the [WHO catalog](https://www.who.int/publications/i/item/9789240028173) but in positions that are different to those of the catalog. It stores the results in a new folder called 'SNPs_in_DRgenes'. This option will not look for SNPs in the catalog, so it will not generate the same outputs as the normal mode of running. 

```bash
ThePipeline2 resistance
```

The only argument for this module is:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -t, --threads | THREADS | Number of threads to use. Default = 1.|
| -ns, --new-snps | | Wether to look only for alternative SNPs in DR genes that are not in the catalog.|

If everything was fine, the resistance module generates the following output:

- resistance_result.tsv: A tabular file with one row per sample. First column correspond to sample ID, whereas the rest of the columns correspond to each antibiotic tested. It contains the resistance profile derived from the comparison of the WHO catalog and the positions called for each sample. The different values that can be found, for each specific drug, are:
  - WT: All positions has been called as wild-type, or if there is a DR allele present in the catalog, is under 0.05 frequency.
  - Unknown: It may happem that:
    - At least one position of the catalog has low coverage.
    - At least one position of the catalog don't present in the PREFIX.snp.vcf, PREFIX.indel.vcf, PREFIX.wt.tsv or PREFIX.lowcov.tsv files.
  - Likely resistant: At least one mutation of the catalog is present in the sample, with a frequency above 0.90.
  - Possible resistant: At least one mutation of the catalog is present in the sample, with a frequency ranging from 0.05 to 0.90.
  - New SNP: One mutation found in the sample is above 0.05 frequency and matches a position having resistant mutations in the catalog, but the allele is not reported in the catalog. For double mutations in the same codon, it will only report new alleles if the codon has been already cataloged as having other double mutations generating resistance (i.e there are other double mutations in that codon in the WHO catalog). 
  - New Indel: One deletion found in the sample is above 0.05 frequency and starts in the same position as other DR deletions in the catalog. However, this exact deletion is not reported in the catalog (has a different size).
- PREFIX.snp_report.tsv: A tabular file only generated for those samples having a profile different to 'WT' or 'Unknown' (one file per sample). In this file, you will find one row per drug and DR mutation. This means that mutations affecting more than one drug will appear more than one time. Fields that can be found are:
  - Drug: Drug to which the mutation potentially confers resistance.
  - Genomic_pos: Genomic coordinate (based on H37Rv) of the mutation. In the case of deletions, it will mark the previous position of the deletion start.
  - Type: Either SNP or Del(etion).
  - WT: Wild-type allele in the MTBC ancestor genome.
  - ALT: Alternative allele found in the .snp.vcf or indel.vcf files.
  - Gene: Gene in which the mutation is found.
  - Freq: Frequency of the alternative allele [0 - 1].
  - N_reads: Number of reads having the alternative allele.
  - Variant_name: Name of the variant in the WHO catalog. It will be empty if the allele is new (not present in the catalog).
  - Codon_change: Nucleotide change provoked in the codon. It will be empty if the allele is new (not present in the catalog).
  - AA_change: Amino acid change provoked in the amin acid sequence. It will be empty if the allele is new (not present in the catalog).
  - Evidence: Evidence supporting this mutation, as it is reported in the WHO catalog. It will report 'New deletion' or 'New allele' if the position is reported in the WHO catalog but he allele found it is not.
  - In the case of new alleles, the codon_change and AA_change will not be written in this file. However, you can find this information in the PREFIX.snp.vcf file, as this file has been annotated in the calling process with SNPEff.
- SNPs_in_DRgenes/PREFIX.SNPs_DRgenes.vcf: Created only when the -ns,--new-snps option is set. This modules creates a new folder called 'SNPs_in_DRgenes'. Inside, for each sample, it writes the rows of the PREFIX.snp.vcf file that correspond to SNPs falling in DR genes but that are not present in the WHO catalog. WARNING: these SNPs may or may NOT be associated with DR! 

**Examples of usage:**

Basic
```bash
ThePipeline2 resistance -t 15
```
Alternative mode of running
```bash
ThePipeline2 resistance -t 15 -ns
```
---

## Lineage typing (typing)
This module compares the different variants identified by the [calling](#variant-calling-calling) module againts a [panel of lineage definitory SNPs](http://tgu.ibv.csic.es/?page_id=1794). It scans all the .snp.vcf files present in the folder in which it is executed. It returns a report with the lineage of each sample. It will also warn of potential mixed infections if any of the lineage definitory SNPs is found in the sample at a frequency under 90%. The panel of SNPs is stored in the /ThePipeline2/data/ folder, and is named [snp_phylo.tsv](http://tgu.ibv.csic.es/wp-content/uploads/2021/05/snps_for_typing.txt). 

```bash
ThePipeline2 typing
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -s, --short| | Use a reduced version of the SNP list. A higher level of classification is obtained. For example: lineage4.1.2 instead of lineage4.1.2.1. By default is not selected|
|-p, --prefix |PREFIX| It will be used for naming the outputs generated. By default current date (YearMonthDay)|

If everything went well, the typing module generates the following output:

- PREFIX.samples_lineages.tsv: A tabular file with one row per sample. First column correspond to sample ID, second column corresponds to the MTBC lineage infered. In the second column you may find:
  - 'Unknown': If no phylogenetic SNP was found in the .snp.vcf files
  - lineageXXX: The deepest phylogenetic MTBC lineage found in the sample. For example, if the script found SNPs definitory of lineage2, lineage2.2 and lineage2.2.4, only lineage2.2.4 will be reported. 
  - 'Possible mixed infection: ': If any of the lineage definitory SNPs is found at a frequency under 90%. After the two points (:), information about the genomic position, the frequency, and the lineage of each phylogenetic definitory SNPs found will be shown.       

**Examples of usage:**

Basic
```bash
ThePipeline2 typing
```

Definig a prefix and using the short version of the SNP panel
```bash
ThePipeline2 typing -s -p run103 
```

---

## Calculate genetic distances (distances)

This module calculates the pairwise genetic distance (in SNP number) between each pair of sequences of a multifasta file. This module needs the [pairsnp-python](https://github.com/gtonkinhill/pairsnp) library to run.

```bash
ThePipeline2 distances
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -f, --fasta| | Multifasta file containing the nucleotide sequences for which the pairwise distances will be calculated|
| -t, --threads | THREADS | Number of threads to use. Default = 1.|
|-p, --prefix |PREFIX| It will be used for naming the outputs generated. By default current date (YearMonthDay)|
|-l, --limit | LIMIT | Generates an extra output, with only the distances that are lower or equal to the value LIMIT specified|

If everything went well, the distances module generates the following output:

- PREFIX.genetic_distances.tsv: A sorted tabular file with one row per sequence. You will find the two sequences of the pair (columns Sequence_1 and Sequence_2) and the genetic distance in SNP number (Distance).
- PREFIX.genetic_distances_LIMIT.tsv: A sorted tabular file with one row per sequence. You will find the two sequences of the pair (columns Sequence_1 and Sequence_2) and the genetic distance in SNP number (Distance). It is generated if the argument _-l, --limit_ is specified. Only distances less or equal to the limit specified are written in this file.

**Examples of usage:**

Basic
```bash
ThePipeline2 distances -t 20 -f 20220726.mf.fasta -p run1
```

---

## Calculate transmission clusters (getclusters)

This module calculates the number of transmission clusters using a defined SNP threshold. It uses as an input the PREFIX.genetic_distances.tsv table generated using [ThePipeline2 distances](#calculate-genetic-distances-distances) module.

```bash
ThePipeline2 getclusters
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -d, --distances| DISTFILE | Genetic distances file. Data should be in the form of a table, either tabular or with comma separated values. Name of the sequences must be in the first and third column, and distance (in SNP number) in the second column. It is prepared to work with the table generated with ThePipeline2 distances. |
| -thres, --threshold | THRESHOLD | SNP threshold for defining recent transmission. By default = 5|
| -sep | [Space\|Tab] | Character separator of the distance file. By default = Tab|
|-osi, --output-senyorito-irving| | Alternative output format. Legacy argument.|
|-p, --prefix |PREFIX| It will be used for naming the outputs generated. By default current date (YearMonthDay)

If everything went well, the getclusters module generates the following output:

- PREFIX.clusters.tsv: By default, a tabular file with one row per transmission cluster, with header. The first column correspond to the Cluster ID, while the second correspond to the samples included in this cluster. If the _-osi_ argument is selected, the table will have one row per sample. The first column correspond to the Cluster ID, while the second column is the Sample ID.

**Examples of usage:**

Basic
```bash
ThePipeline2 getclusters -d 20220907.genetic_distances.tsv
```


---

## Generate MultiQC report (multiqc)
This module uses [MultiQC](https://multiqc.info/) to generate a general report, joining partial reports derived from some of the tools used in the other modules execution. It joins the reports derived from [fastclean](#fastq-filtering-fastclean), [kraken](#taxonomic-filtering-kraken) and [mapping](#mapping-mapping), creating an interactive HTML webpage with information about the trimming process performed with fastp, the contaminations detected with kraken, the potential PCR and optical duplicates detected with PICARD tools.
```bash
ThePipeline2 multiqc
```

The main arguments for this module are:

|Argument|Value|Description|
|:------------- |:---|:-----|
| -o, --output| OUTPUT | Output file name. By default will be multiqc_report. If this file alreay exists in the folder, it will be overwritten |
| -f, --folder |FOLDER| Folder containing the logs and reports that MultiQC will use as input. If you have previously run _ThePipeline2 organize_, they will be in the reports/ folder|

If everything was fine, the multiqc module generates the following output:
- OUTPUT.html: An interactive HTML file with information about the trimming process performed with fastp, the contaminations detected with kraken and the potential PCR and optical duplicates detected with PICARD tools

- OUTPUT_data/: A folder containing all the information needed by OUTPUT.html to be properly visualized in a web browser. This folder needs to be in the same location as OUTPUT.html.

**Examples of usage:**

Basic
```bash
ThePipeline2 multiqc -o 220401run_multiqc
```
---

