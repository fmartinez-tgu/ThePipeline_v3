#! /usr/bin/env python3.7
# Copyrigth (C) 2022 Alvaro Chiner Oms

# This file is part of ThePipeline3
# ThePipeline3 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# ThePipeline3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with ThePipeline3.  If not, see <http://www.gnu.org/licenses/>.


def Organize():
    import os
    import subprocess as sp

    # list of folder to be created
    folders = ["reports", "original_fastq", "clean_fastq",
               "filtered_fastq", "BAM_CRAM", "variants",
               "resistance"]

    for folder in folders:
        # check if any folder exists
        # and ask for confirmation
        # Files could be overwritten!
        if os.path.isdir(folder):
            ans = ""
            while ans != 'y' and ans != 'n':
                ans = input("\033[93mWARNING: Some folders already exist. "
                            "Data will be moved inside existing folders, "
                            "and files with same names will be overwritten. "
                            "Proceed? (y/n):\033[0m")
            if ans == 'n':
                os._exit(0)  # if 'n' stop execution
            elif ans == 'y':
                break  # if 'y', break loop and continue

    for folder in folders:
        try:
            os.mkdir(folder)  # create folders
        except OSError:
            print(folder, "\033[93mWARNING: folder already exists\033[0m")

    # Move files to reports/
    sp.run("mv *.kraken.report reports/", shell=True, capture_output=True)
    sp.run("mv *.cleanlog reports/", shell=True, capture_output=True)
    sp.run("mv *.coverage reports/", shell=True, capture_output=True)
    sp.run("mv *.dup.metrix reports/", shell=True, capture_output=True)
    sp.run("mv *.json reports/", shell=True, capture_output=True)
    sp.run("mv *.contaminants reports/", shell=True, capture_output=True)
    sp.run("mv *.meancov reports/", shell=True, capture_output=True)
    sp.run("mv *multiqc* reports/", shell=True, capture_output=True)
    sp.run("mv *.history reports/", shell=True, capture_output=True)
    sp.run("mv *_stats reports/", shell=True, capture_output=True)

    # Move files to BAM_CRAM/
    sp.run("mv *.bam BAM_CRAM/", shell=True, capture_output=True)
    sp.run("mv *.cram BAM_CRAM/", shell=True, capture_output=True)

    # Move fastq files
    sp.run("mv *1.fastq.gz original_fastq/", shell=True,
           capture_output=True)
    sp.run("mv *2.fastq.gz original_fastq/", shell=True,
           capture_output=True)
    sp.run("mv *P1.clean.fastq.gz clean_fastq/", shell=True,
           capture_output=True)
    sp.run("mv *P2.clean.fastq.gz clean_fastq/", shell=True,
           capture_output=True)
    sp.run("mv *P1.filtered.fastq.gz filtered_fastq/", shell=True,
           capture_output=True)
    sp.run("mv *P2.filtered.fastq.gz filtered_fastq/", shell=True,
           capture_output=True)

    # Move files to variants/
    sp.run("mv *.vcf variants/", shell=True, capture_output=True)
    sp.run("mv *.annoF variants/", shell=True, capture_output=True)
    sp.run("mv *.inferred_WT.tsv variants/", shell=True, capture_output=True)
    sp.run("mv *.lowcov.tsv variants/", shell=True, capture_output=True)

    # Move files from the resistance analysis
    sp.run("mv *.snp_report.tsv resistance/", shell=True, capture_output=True)
    sp.run("mv *resistance_result.tsv resistance/",
           shell=True, capture_output=True)
