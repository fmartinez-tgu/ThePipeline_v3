#! /usr/bin/env python3.7
# Copyright (C) 2025 Miguel Moreno Molina & Francisco Jose Martinez

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

# Module for variant calling, using Varscan2, Mutect2 and Minos

import vcf
def check_extension(ext):
    '''Check if the extension is either BAM or CRAM'''
    import sys

    if ext.endswith(".bam") or ext.endswith(".cram"):
        return 0
    else:
        sys.exit("\033[91mERROR: We only accept files"
                 " ended in .bam or .cram extension\033[0m")

def get_reference_id(prefix, ext, samtools):
    '''Function to extract the name of the reference to be added to the output SNP files, so it's not MTB_anc by default'''
    import subprocess
    ref_ID = subprocess.run([f"{samtools} view -H {prefix}{ext} | grep '^@SQ' | cut -f 2 | cut -f 2 -d ':'"], stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    ref_ID = ref_ID.stdout.strip().split('\n')
    ref_name = subprocess.run([f"{samtools} view -H {prefix}{ext} | grep '^@SQ' | cut -f 5 | rev | cut -f 1 -d '/' | rev"], stdout=subprocess.PIPE, shell=True, universal_newlines=True)
    ref_name = ref_name.stdout.strip()
    ref_ID.append(ref_name)
    return ref_ID

def VCFtoPandas(file):
    import pandas
    import io

    with open(file, 'r') as f:
        lines = [li for li in f if not li.startswith('##')]
    data = pandas.read_csv(io.StringIO(''.join(lines)),
                           dtype={'#CHROM': str, 'POS': int, 'ID': str,
                                  'REF': str, 'ALT': str, 'QUAL': str,
                                  'FILTER': str, 'INFO': str},
                           sep='\t')

    return data


def mapq60_filter(prefix, ext, samtools, reference):
    '''Before running VarScan, we're going to filter the CRAM file using a MAPQ >= 60 to use it exclusively with VarScan'''
    print(f"FILTERING MAPQ 60 IN INPUT BAM")
    import subprocess as sp
    from subprocess import call

    if ext == ".sort.bam":
        bamfile = f"{prefix}{ext}"
        # Command to filter the BAM file using a mapq >= 60 
        cmd = f"""{samtools} view -h {bamfile} | awk '{{ if ($1 ~ /^@/ || $5 >= 60) print }}' > {prefix}_mapq60.sort.bam"""
        sp.run(cmd, stdout=sp.PIPE, shell=True, universal_newlines=True)

    elif ext == ".cram":
        cram_file = f"{prefix}{ext}"
        cram_to_sortbam = [samtools, "view", "-T", reference, "-B", "-o", "{}.sort.bam".format(prefix), cram_file]
        call(cram_to_sortbam)

        bamfile = f"{prefix}.sort.bam"
        cmd = f"""{samtools} view -h {bamfile} | awk '{{ if ($1 ~ /^@/ || $5 >= 60) print }}' > {prefix}_mapq60.sort.bam"""
        sp.run(cmd, stdout=sp.PIPE, shell=True, universal_newlines=True)


def VarScan(reference, prefix, varscan, samtools, ext, ref_ID):
    '''Call samtools mpileup and Varscan, then transform output to simple VCF for Minos'''
    from subprocess import call
    import sys
    import pandas as pd
    from .History import UpdateHistory

    mapq60_filter(prefix, ext, samtools, reference)

    bam_filtered_file = f"{prefix}_mapq60.sort.bam"

    with open("{}.mpileup.remove".format(prefix), "w") as outfh:
        cmd = [samtools, "mpileup", "-AB", "-f", reference, bam_filtered_file]
        stat = call(cmd, stdout=outfh)
        if stat != 0:
            sys.exit("Pipeline stopped at samtools mpileup!\n")
        UpdateHistory(cmd, "samtools", prefix)


    # VarScan calling
    with open("{}.snp".format(prefix), "w") as outfh:
        cmd = ["java", "-Xms10G", "-Xmx32G", "-jar", varscan,
            "pileup2snp", "{}.mpileup.remove".format(prefix),
            "--min-coverage", "3",
             "--min-reads2", "3",
             "--min-freq-for-hom", "0.90",
             "--min-var-freq", "0.05"]
        stat = call(cmd, stdout=outfh)
        if stat != 0:
            sys.exit("Pipeline stopped at VarScan pileup2snp!\n")
        UpdateHistory(cmd, "varscan", prefix)
        call(["rm", "{}.mpileup.remove".format(prefix)])

        UpdateHistory(cmd, "varscan", prefix)

    # Conversion of VarScan output to VCF so it can be parsed by Minos
    with open("{}.parsed.vcf".format(prefix), "w") as outfh:
        outfh.write("##fileformat=VCFv4.2\n")
        outfh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        outfh.write(f"##contig=<ID={'_'.join(ref_ID[0:-1])},length=4411532,assembly={ref_ID[-1]}>\n")
        outfh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + prefix + "\n")
        snpfile = pd.read_csv("{}.snp".format(prefix), sep='\t', header=0)
        for index, snp in snpfile.iterrows():
            depth = str(int(snp["Reads1"]) + int(snp["Reads2"]))
            f = str(float(snp["VarFreq"].replace("%","").replace(",","."))/100)
            outfh.write("\t".join([snp["Chrom"], str(snp["Position"]), ".", snp["Ref"],
                        snp["VarAllele"], ".", ".", ".", "GT:AD:AF", "0/1:" + depth + ":" + f + "\n"]))


def Mutect2(reference, prefix, gatk, samtools, genomeCoverageBed,
            ext, threads, min_depth, min_qual, freq):
    '''Call Mutect2 to perfom variant calling and obtain a VCF file
    and bedtools to calculate coverage and determine lowcov positions'''
    import subprocess as sp
    import sys
    import os
    import pandas
    from .History import UpdateHistory

    # need to import libstdc++.so.6 updated library for
    # bedtools to work. It's in /PipeModules/Configs/lib
    dirname = os.path.split(os.path.dirname(os.path.realpath(__file__)))[0]
    os.environ['LD_LIBRARY_PATH'] = os.path.join(dirname,
                                                 'data',
                                                 'libs')
    

    # index BAM
    cmd_indexBam = [samtools, "index", "-@", threads,
                    "{}.sort.bam".format(prefix)]
    stat = sp.call(cmd_indexBam)
    UpdateHistory(cmd_indexBam, "samtools", prefix)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: "
                 "Pipeline stopped when indexing BAM!\n\033[0m".format(prefix))

    # CALLING WITH MUTECT2
    # Run Mutect2
    cmd_Mutect = [gatk, "Mutect2",
                   "-R", reference, "-I", "{}.sort.bam".format(prefix),
                   "-O", "{}_unfiltered.vcf".format(prefix), "-OVI", "false",
                   "--verbosity", "ERROR",
                   "--QUIET", "true", "-mbq", min_qual,
                   "--callable-depth", min_depth,
                   "-mnp-dist", "0",  # important to not filter MNPs later,
                   "--linked-de-bruijn-graph", "true",
                   "--native-pair-hmm-threads", threads, "--f1r2-tar-gz", f"{prefix}.f1r2.tar.gz"]
    
    stat = sp.call(cmd_Mutect)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: Pipeline stopped"
                 " when running Mutect2!\n\033[0m".format(prefix))
    UpdateHistory(cmd_Mutect, "gatk", prefix)
    
    # Now we mark the reads with orientation bias and filter them
    print("\033[92m\nFiltering {} by read orientation\n\033[00m".format(prefix))
    sp.run(f"{gatk} LearnReadOrientationModel -I {prefix}.f1r2.tar.gz -O read-orientation-model_{prefix}.tar.gz", shell=True)
    sp.run(f"{gatk} FilterMutectCalls -V {prefix}_unfiltered.vcf -R {reference} --ob-priors read-orientation-model_{prefix}.tar.gz -O {prefix}.vcf --min-median-read-position 15", shell=True)

    with open(f"{prefix}.vcf", "r+") as input_file:
        lines = input_file.readlines()
        headers = [x.strip() for x in lines if '#' in x]
        lines_noheaders = [x.strip() for x in lines if "#" not in x]

    with open(f"{prefix}_no_orientation.vcf", "w+") as output_file:
        output_file.write("\n".join(headers))

        for line in lines_noheaders:
            if "orientation" not in line:
                output_file.write(f"\n{line}")


    #now produce the gVCF for accurate WT calling
    print("\033[92m\nObtaining {}.gvcf\n\033[00m".format(prefix))
    cmd_Mutect2 = [gatk, "Mutect2",
                   "-R", reference, "-I", "{}.sort.bam".format(prefix),
                   "-O", "{}.gvcf".format(prefix), "-OVI", "false",
                   "--verbosity", "ERROR",
                   "--QUIET", "true", "-mbq", min_qual,
                   "--callable-depth", min_depth,
                   "-mnp-dist", "0",  # important to not filter MNPs later
                   "--linked-de-bruijn-graph", "true",
                   "--native-pair-hmm-threads", threads,
                   "-ERC", "BP_RESOLUTION"]
    
    stat = sp.call(cmd_Mutect2)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: Pipeline stopped"
                 " when running Mutect2 for gVCF generation!\n\033[0m".format(prefix))
    UpdateHistory(cmd_Mutect2, "gatk", prefix)

    # Run SelectVariants to split INDELs and SNPs
    cmd_selectI = [gatk, "SelectVariants", "-R", reference,
                   "-V", "{}_no_orientation.vcf".format(prefix),
                   "--verbosity", "ERROR",
                   "--QUIET", "true",
                   "-OVI", "false",
                   "-select-type", "INDEL", "-O",
                   "{}.indel.vcf".format(prefix)]
    stat = sp.call(cmd_selectI)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: "
                 "Pipeline stopped at SelectVariants\n\033[0m".format(prefix))
    UpdateHistory(cmd_selectI, "gatk", prefix)

    cmd_selectS = [gatk, "SelectVariants", "-R", reference,
                   "-V", "{}_no_orientation.vcf".format(prefix),
                   "--verbosity", "ERROR",
                   "--QUIET", "true",
                   "-OVI", "false",
                   "-xl-select-type", "INDEL",
                   "-xl-select-type", "MNP",
                   "-xl-select-type", "SYMBOLIC",
                   "-xl-select-type", "NO_VARIATION",
                   "-O", "{}.snp.vcf".format(prefix)]
    stat = sp.call(cmd_selectS)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: "
                 "Pipeline stopped at SelectVariants\n\033[0m".format(prefix))
    UpdateHistory(cmd_selectS, "gatk", prefix)

    # Extract multiallelic from indel.vcf
    sp.run("grep -v 0/1/2 {0}.indel.vcf > {0}.indel_sin_multiallelic.vcf".format(prefix), shell=True)
    UpdateHistory("grep -v 0/1/2 {0}.indel.vcf > {0}.indel_sin_multiallelic.vcf".format(prefix), "custom", prefix)

    # Filter indels 
    indel_reader = vcf.Reader(open("{}.indel_sin_multiallelic.vcf".format(prefix), "r"))
    indel_vcf_remade = vcf.Writer(open("{}.remade.indel.vcf".format(prefix), "w+"), indel_reader)

    for record in indel_reader:
        if(record.samples[0]['AF'] >= 0.05 and record.samples[0]['DP'] >= 10):
            indel_vcf_remade.write_record(record)

    indel_vcf_remade.close()

    # CALCULATE LOWCOV FROM COVERAGE
    # Calculate coverage
    with open("{}.coverage".format(prefix), "w") as outfh:
        cmd_cov = [genomeCoverageBed, "-ibam", "{}.sort.bam".format(prefix),
                   "-d", "-g", reference]
        stat = sp.call(cmd_cov, stdout=outfh)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: Pipeline stopped when"
                 " calculating BAM coverage!\n\033[0m".format(prefix))
    outfh.close()
    UpdateHistory(cmd_cov, "bedtools", prefix)

    # lowcov will be all bases bellow min_depth cut-off
    cov_info = pandas.read_csv("{}.coverage".format(prefix),
                               sep='\t')
    cov_info.columns = ['CHROM', 'POS', 'COV']
    cov_info[cov_info['COV'] < float(min_depth)
             ]['POS'].to_csv("{}.lowcov.tsv".format(prefix),
                             sep='\n',
                             index=False)


# Now we filter VarScan (both parsed and original) and Mutect2 outputs using the annotation filter

def annoF_varscan_mutect(prefix):

    import os
    annotation = LoadAnnotation("/data/ThePipeline_v3/data/H37Rv.annotation_new.tsv".format(
                    os.path.split(
                        os.path.dirname(
                            os.path.abspath(__file__)))[0]))

    # VarScan files
    with open("{}.parsed.vcf_annoF".format(prefix), "w") as outfile:
        with open("{}.parsed.vcf".format(prefix)) as infile:
            lines_varscan = infile.readlines()
            headers = [x.strip() for x in lines_varscan if "#" in x]
            rest_of_file = [x.strip() for x in lines_varscan if "#" not in x]
            outfile.write("\n".join(headers)+"\n")
            for line in rest_of_file:
                pos = line.strip().split("\t")[1]
                if pos not in annotation:
                        outfile.write(f"{line}\n")

    with open("{}.snp_annoF".format(prefix), "w") as outfile:
        with open("{}.snp".format(prefix)) as infile:
            header = infile.readline()
            outfile.write(header)
            for line in infile:
                pos = line.strip().split("\t")[1]
                if pos not in annotation:
                        outfile.write(line)

    # Mutect2 file
    with open("{}.snp.vcf_annoF".format(prefix), "w") as outfile:
        with open("{}.snp.vcf".format(prefix)) as infile:
            lines_mutect2 = infile.readlines()
            headers = [x.strip() for x in lines_mutect2 if '#' in x]
            rest_of_file = [x.strip() for x in lines_mutect2 if '#' not in x]
            outfile.write("\n".join(headers)+"\n")
            for line in rest_of_file:
                pos = line.strip().split("\t")[1]
                if pos not in annotation:
                    outfile.write(f"{line}\n")

    # Remove the original file and just keep the one filtered by annotation
    os.rename(f"{prefix}.parsed.vcf", f"{prefix}.parsed.vcf.original_no_annoF")
    os.rename(f"{prefix}.parsed.vcf_annoF", f"{prefix}.parsed.vcf")
    os.rename(f"{prefix}.snp", f"{prefix}.snp.original_no_annoF")
    os.rename(f"{prefix}.snp_annoF", f"{prefix}.snp")
    os.rename(f"{prefix}.snp.vcf", f"{prefix}.snp.vcf.original_no_annoF")
    os.rename(f"{prefix}.snp.vcf_annoF", f"{prefix}.snp.vcf")


def mutect2_vcf_to_tab(prefix):

    '''Convert Mutect2 VCF output to tab file for easier visualization. Original VCF will be removed unless -kmvcf parameter is used'''
    # Import Mutect2 VCF file
    
    with open(f"{prefix}.snp.mutect", "r+") as mutect_file:
        lines_mutect = mutect_file.readlines()

    lines_mutect = [x for x in lines_mutect if '#' not in x] # Remove all headers

    # Start up new file
    with open(f"{prefix}.snp.mutect.tab", "w+") as output_mutect_tab:
        output_mutect_tab.write("#Chrom\tPosition\tRef\tCons\tVarFreq\tCov_allele\tVarAllele\n")
    
        # Now, we convert the Mutect2 lines to the new format and write them in the output file
        for line in lines_mutect:
            tokens = line.split("\t")
            varfreq = str(round(float(tokens[9].split(":")[2])*100,2)) # Save the variant frequency 
            cov_allele = tokens[9].split(":")[3] # Save the depth 
            new_line = [tokens[0],tokens[1],tokens[3],tokens[4],varfreq,cov_allele,tokens[4]] # Line that will be added to the file
            output_mutect_tab.write(("\t").join(new_line)+"\n")



# Filter multiallelic from Mutect2 snp output
def filter_multiallelic_from_mutect2_snp(prefix):
    from .History import UpdateHistory
    '''Filter multiallelic positions from Mutect2 SNP VCF output, saving them in a separate file and remaking the original VCF 
       so that each multiallelic SNP is split into multiple rows. Multiallelic positions will be saved as two different SNPs in the remade VCF 
       so that Minos can take them into account.'''

    with open("{}.snp.vcf".format(prefix), "r+") as input_file:
        lines = input_file.readlines()
    
    headers = [line for line in lines if line.startswith("#")]
    lines_header_CHROM = [index for index, element in enumerate(headers) if "#CHROM" in element]
    lines_sin_headers = lines[int(lines_header_CHROM[0])+1:]

    # File to save multiallelic positions
    multiallelic_vcf = open("{}.multiallelic.snp.vcf".format(prefix), "w+")
    multiallelic_positions = []


    with open("{}.remade.snp.vcf".format(prefix), "w+") as output_file:
        output_file.write("".join(headers))

        for line in lines_sin_headers: # For each position
            tokens = line.split("\t")
            
            ALT = tokens[4]
            if "," not in ALT: # If it's not multiallelic, we write it directly to the output file
                output_file.write(line)
            else:
                multiallelic_positions.append(line.strip())

                # First, we check if any of the possible variants of the position is not a SNP
                # We get the length of the possible variants and if there's any non-SNP variants (chunks) we discard it. Since the indels are saved in another file

                ALT_list = ALT.split(",")
                lengths = [len(element) for element in ALT_list] 
                lengths_greater_than_1 = [i for i in lengths if i > 1] 
                
                if len(lengths_greater_than_1) > 0: 
                    continue
                
                else: 
                    
                    AS_FilterStatus_to_write = tokens[7].split(";")[0]
                    AS_SB_TABLE_info_list = tokens[7].split(";")[1][12:].split("|")
                    AS_SB_TABLE_info_list_list_to_write = []
                    for i in range(1,len(AS_SB_TABLE_info_list)):
                        AS_SB_TABLE_info_list_list_to_write.append([AS_SB_TABLE_info_list[0],AS_SB_TABLE_info_list[i]])
                    
                    DP_to_write = tokens[7].split(";")[2]
                    ECNT_to_write = tokens[7].split(";")[3]

                    GERMQ_to_write = tokens[7].split(";")[4]

                    MBQ_list = tokens[7].split(";")[5][4:].split(",")
                    MBQ_list_to_write = []
                    for i in range(1,len(MBQ_list)):
                        MBQ_list_to_write.append([MBQ_list[0],MBQ_list[i]])
                    
                    MFRL_list = tokens[7].split(";")[6][5:].split(",")
                    MFRL_list_to_write = []
                    for i in range(1,len(MFRL_list)):
                        MFRL_list_to_write.append([MFRL_list[0],MFRL_list[i]])
                    
                    MMQ_list = tokens[7].split(";")[7][4:].split(",")
                    MMQ_list_to_write = []
                    for i in range(1,len(MMQ_list)):
                        MMQ_list_to_write.append([MMQ_list[0],MMQ_list[i]])

                    MPOS_list_to_write = tokens[7].split(";")[8][5:].split(",")
                    POPAF_list_to_write = tokens[7].split(";")[9][6:].split(",")
                    ROQ_to_write = tokens[7].split(";")[10]
                    TLOD_list_to_write = tokens[7].split(";")[11][5:].split(",")

                    GT_to_write = "0/1"

                    AD_list = tokens[9].split(":")[1].split(",")
                    AD_list_to_write = []
                    for i in range(1,len(AD_list)):
                        AD_list_to_write.append([AD_list[0],AD_list[i]])
                    
                    AF_list_to_write = tokens[9].split(":")[2].split(",")
                    DP_to_write = tokens[9].split(":")[3]

                    F1R2_list = tokens[9].split(":")[4].split(",")
                    F1R2_list_to_write = []
                    for i in range(1,len(F1R2_list)):
                        F1R2_list_to_write.append([F1R2_list[0],F1R2_list[i]])

                    F2R1_list = tokens[9].split(":")[5].split(",")
                    F2R1_list_to_write = []
                    for i in range(1,len(F2R1_list)):
                        F2R1_list_to_write.append([F2R1_list[0],F2R1_list[i]])

                    FAD_list = tokens[9].split(":")[6].split(",")
                    FAD_list_to_write = []
                    for i in range(1, len(FAD_list)):
                        FAD_list_to_write.append([FAD_list[0],FAD_list[i]])
                    
                    SB_to_write = tokens[9].split(":")[7]

                    for i in range(0,len(ALT_list)):
                        output_file.write('''{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7};AS_SB_TABLE={8};{9};{10};{11};MBQ={12};MFRL={13};MMQ={14};MPOS={15};POPAF={16};{17};TLOD={18}\t{19}\t{20}:{21}:{22}:{23}:{24}:{25}:{26}:{27}'''.format(
                            tokens[0],
                            tokens[1],
                            tokens[2],
                            tokens[3],
                            ALT_list[i],
                            tokens[5],
                            tokens[6],
                            AS_FilterStatus_to_write,
                            ("|").join(AS_SB_TABLE_info_list_list_to_write[i]), 
                            DP_to_write,
                            ECNT_to_write,
                            GERMQ_to_write,
                            (",").join(MBQ_list_to_write[i]),
                            (",").join(MFRL_list_to_write[i]),
                            (",").join(MMQ_list_to_write[i]),
                            MPOS_list_to_write[i],
                            POPAF_list_to_write[i],
                            ROQ_to_write,
                            TLOD_list_to_write[i], 
                            tokens[8], 
                            GT_to_write, 
                            (",").join(AD_list_to_write[i]),
                            AF_list_to_write[i],
                            DP_to_write,
                            (",").join(F1R2_list_to_write[i]),
                            (",").join(F2R1_list_to_write[i]),
                            (",").join(FAD_list_to_write[i]),
                            SB_to_write
                        ))           
    
    # We save the multiallelic positions in a separate VCF file
    multiallelic_vcf.write("".join(headers))
    multiallelic_vcf.write("\n".join(multiallelic_positions))
    multiallelic_vcf.close()
    UpdateHistory("Multiallelic filtering and trialelic readjudication for Mutect2 output", "custom", prefix)
    return 0


def callWT(prefix):
    '''Creates a Wild Type (WT) position file along with their read depth'''
    import vcf
    from .History import UpdateHistory
    vcf_reader = vcf.Reader(open("{}.gvcf".format(prefix), 'r'))

    wt = {}
    for record in vcf_reader:
        position, ref, alt  = record.POS, record.REF, record.ALT
        depth, genotype = record.samples[0]["DP"], record.samples[0]["GT"]
        alleles_depth = record.samples[0]["AD"]

        if genotype == "0/0" and float(alleles_depth[0]) > (float(depth * 0.90)):
            wt[position] = depth

	#if genotype is not 0/0, again make sure there are no variants above 10%
	#frequency, and then report the correct depth for the WT allele

        elif len(alt) > 1:
            frequencies = record.samples[0]["AF"]

            if all([True if freq <= 0.10 else False for freq in frequencies]):
                wt[position] = alleles_depth[0]

    # save WT positions
    with open("{}.wt".format(prefix), 'w') as WT_file:
        WT_file.write("Position\tWT_Depth\n")
        for pos, d in wt.items():
            if d >= 3:
                WT_file.write(str(pos) + "\t" + str(d) + "\n")
    WT_file.close()
    UpdateHistory("WT obtained", "custom", prefix)



def Minos(reference,minos, snpEff, prefix, single_end, ref_ID):
    '''Perform final variant call adjudication with Minos'''
    import sys, os
    import subprocess as sp
    from subprocess import call
    from .History import UpdateHistory
    
    if not single_end:
        cmd = ["singularity", "exec", "--bind", "/"+os.getcwd().split("/")[1],
            minos, "minos", "adjudicate", "--force",
            "--reads", "{}.P1.filtered.fastq.gz".format(prefix),
            "--reads", "{}.P2.filtered.fastq.gz".format(prefix),
            "{}_minos".format(prefix), reference.split("/")[-1],
            "{}.remade.snp.vcf".format(prefix), "{}.parsed.vcf".format(prefix)]
    else:
         cmd = ["singularity", "exec", "--bind", "/"+os.getcwd().split("/")[1],
            minos, "minos", "adjudicate", "--force",
            "--reads", "{}.filtered.fastq.gz".format(prefix),
            "{}_minos".format(prefix), reference.split("/")[-1],
            "{}.remade.snp.vcf".format(prefix), "{}.parsed.vcf".format(prefix)]       
           
    stat = call(cmd)
    if stat != 0:
        sys.exit("Pipeline stopped at Minos variant adjudication for sample {}!\n".format(prefix))
    UpdateHistory(cmd, "minos", prefix)

    # We extract the Minos output from each folder and filter WT positions
    sp.run(f"grep -v '0/0' {prefix}_minos/final.vcf > {prefix}.final_sin_wt.vcf; rm -r {prefix}_minos", shell=True)
    UpdateHistory("WT filtered from Minos output", "custom", prefix)


def minos_raw_vcf_to_tab(prefix, ref_ID, snpEff):
    """
    Convert Minos consensus VCF and complementary caller outputs into a tabular summary.
    Summary:
        Reads Minos consensus VCF and the outputs from VarScan and Mutect2 (expected to follow specific filename patterns
        derived from the provided prefix). It complements Minos consensused positions with positions that are common between
        VarScan and Mutect2 but missing from Minos, then builds a tab-delimited summary table that contains, for each
        position: chromosome, position, reference base, consensus call, variant frequency(s) (percent), coverage/depth(s)
        and reported variant allele(s). The function infers frequencies and depths by collecting values from VarScan and
        Mutect2 outputs and averaging when both callers report a value for the same variant. It applies simple heuristics
        to decide how to represent the consensus (IUPAC code, explicit allele, or '?') based on variant frequencies and
        multiplicity. In short:

        - Reads Minos VCF and identifies positions called by VarScan and Mutect2 but missing from Minos.
        - Merges these missing positions into the Minos VCF data.
        - For each variant, collects per-allele frequencies and depths from VarScan and Mutect2,
          taking the mean when both callers report the same allele.
        - Writes a tab file "{prefix}.minos.raw.tab" with columns: Chrom, Position, Ref, Cons, VarFreq, Cov_allele, VarAllele
        - For biallelic sites with ALT freq < 90%: map Ref+ALT to IUPAC code in Cons and keep ALT in VarAllele.
        - For biallelic sites with ALT freq >= 90%: put the ALT in Cons and VarAllele (fixed).
        - For multi-allelic sites:
          - If any allele has >= 90% frequency, that allele is placed in Cons and all alleles/frequencies/depths
            are recorded in VarAllele/VarFreq/Cov_allele fields (major allele first).
          - If no allele is fixed, Cons is set to "?" and alleles/frequencies/depths are written in parallel order
            so frequencies map to the listed variants.
          - Calls filter_raw_minos(...) at the end to apply additional filters (e.g., low coverage or low frequency).
    """

     
    import subprocess as sp
    import pandas as pd
    from subprocess import call

    iupac = {"M": ["A", "C"],
        "R": ["A", "G"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "Y": ["C", "T"],
        "K": ["G", "T"]}
    
    def find_key_with_elements(d, elements):
        '''Finds a key in dictionary d whose value contains all elements in the provided list.'''
        required_elements = set(elements)
        for key, value in d.items():
            if required_elements.issubset(set(value)):
                return key
        return None

    def get_freqs_from_varscan_and_mutect(pos_minos, ref_ID):
        '''This function will obtain the frequency for each consensus position called by Minos from the VarScan and Mutect2 outputs. 
        If the position is present in both callers, the mean frequency will be calculated'''
        
        variantes = cons.split(",") # If there are multiple consensus variants, split them. For each one we will get the frequency from VarScan and Mutect2
        lista_keys = [ref+i for i in variantes]
        dic_variants = {key: [] for key in lista_keys}

        # Get the frequencies from VarScan and Mutect2
        with open("{}.parsed.vcf".format(prefix), "r+") as input_varscan:
            lines_varscan = input_varscan.readlines()
            lines_varscan = [x for x in lines_varscan if "#" not in x]

            for line_varscan in lines_varscan:

                try:
                    if any(f"{item}\t{pos_minos}\t" in line_varscan for item in ref_ID): # If the position is present, add the frequency to the corresponding dictionary key
                        ref_varscan, cons_varscan, freq_varscan = itemgetter(3,4,9)(line_varscan.split("\t"))
                        freq_variant_varscan = freq_varscan.split(":")[-1] # Save the frequency of the variant
                        dic_variants[ref_varscan+cons_varscan].append(float(freq_variant_varscan)) # Save variant frequency in its key
                except KeyError: # Some positions may have 2 possible variants in VarScan. 
                                 # However, if any of them is not present in the ALT column of Minos, it will not be considered
                    continue

                    # Same for Mutect2
        with open("{}.remade.snp.vcf".format(prefix), "r+") as input_mutect2:
            lines_mutect2 = input_mutect2.readlines()

            for line_mutect2 in lines_mutect2:
                try:
                    if any(f"{item}\t{pos_minos}\t" in line_mutect2 for item in ref_ID):
                        ref_mutect, cons_mutect, info_mutect = itemgetter(3,4,9)(line_mutect2.split("\t"))
                        freq_variant_mutect = info_mutect.split(":")[2] 
                        dic_variants[ref_mutect+cons_mutect].append(float(freq_variant_mutect))
                except KeyError:
                    continue
                    
        from statistics import mean

        # Save the mean of the frequencies of each variant extracted from VarScan and Mutect2 and add them 
        # to the VarFreq column ordered from highest to lowest

        freqs_to_write = [round(mean(values)*100,2) for keys,values in dic_variants.items()]
        if len(freqs_to_write) == 1:
            return freqs_to_write[0]
                    
        elif len(freqs_to_write) > 1:
            return freqs_to_write


    def get_depths_from_varscan_and_mutect(pos_minos, ref_ID):
        '''This function will obtain the depth for each consensus position called by Minos from the VarScan and Mutect2 outputs. 
        If the position is present in both callers, the mean depth will be calculated.'''
        
        variantes = cons.split(",")
        lista_keys = [ref+i for i in variantes]
        dic_variants = {key: [] for key in lista_keys}

        # Get depths from VarScan and Mutect2
        with open("{}.parsed.vcf".format(prefix), "r+") as input_varscan:
            lines_varscan = input_varscan.readlines()
    
            for line_varscan in lines_varscan:

                try:
                    if any(f"{item}\t{pos_minos}\t" in line_varscan for item in ref_ID): 
                        ref_varscan, cons_varscan, freq_varscan = itemgetter(3,4,9)(line_varscan.split("\t"))
                        depth_variant_varscan = freq_varscan.split(":")[1] 
                        dic_variants[ref_varscan+cons_varscan].append(float(depth_variant_varscan)) 
                except KeyError: 
                    continue

        with open("{}.remade.snp.vcf".format(prefix), "r+") as input_mutect2:
            lines_mutect2 = input_mutect2.readlines()

            for line_mutect2 in lines_mutect2:
                try:
                    if any(f"{item}\t{pos_minos}\t" in line_mutect2 for item in ref_ID):
                        ref_mutect, cons_mutect, info_mutect = itemgetter(3,4,9)(line_mutect2.split("\t"))
                        depth_variant_mutect = str(int(info_mutect.split(":")[1].split(",")[1])+int(info_mutect.split(":")[1].split(",")[0])) # Guardamos la depth de ambas variantes sumadas
                        dic_variants[ref_mutect+cons_mutect].append(float(depth_variant_mutect))
                except KeyError:
                    continue
                    
        from statistics import mean

        depths_to_write = [int(mean(values)) for keys,values in dic_variants.items()]

        if len(depths_to_write) == 1:
            return depths_to_write[0]
                    
        elif len(depths_to_write) > 1:
            return depths_to_write


    def check_bigger_than_90(list1, val):

        for x in list1:
            if x >= val:
                return True
    
    def check_lower_than_90(list1, val):
     
    # traverse in the list
        for x in list1:
            if x < val:
                return True
    
    from operator import itemgetter

    ### Now we're going to complement the positions consensed by Minos with the common variants between VarScan_Mutect2.

    # Import the original file obtained from Minos with the consensual variants
    with open("{}.final_sin_wt.vcf".format(prefix), "r+") as filtered_raw:
        lines = filtered_raw.readlines()
        lines_headers = [i.strip() for i in lines if "#" in i]
        header = lines_headers[-1]
        lines_preheader = lines_headers[:-1]
        lines_noheader = [i.strip() for i in lines if "#" not in i]

    # We first import the varscan, mutect and minos positions to include the missing ones

    varscan_file = pd.read_csv("{}.snp".format(prefix), sep = "\t", header=0)
    mutect_file = VCFtoPandas("{}.remade.snp.vcf".format(prefix)) # Now .remade.snp.vcf, before was .snp.vcf
    minos_file = VCFtoPandas("{}.final_sin_wt.vcf".format(prefix))
    varscan_dict = varscan_file.set_index('Position')['Cons'].to_dict()
    mutect_dict = mutect_file.set_index('POS')['ALT'].to_dict()
    minos_dict = minos_file.set_index('POS')['ALT'].to_dict()

    # Now, we must obtain all the unique list of positions combining varscan, minos and mutect

    # Get the keys from each dictionary as sets
    varscan_positions = set(varscan_dict.keys())
    mutect_positions = set(mutect_dict.keys())
    minos_positions = set(minos_dict.keys())

    # Common positions
    common_varscan_mutect = varscan_positions & mutect_positions

    # We add the common positions to the final file lines. We only need the MTB_anc, position, ref and alt columns. The rest will be empty
        
    for position in common_varscan_mutect: 

        if position not in minos_positions:
            # Check if the position has only one variant in both callers
            variants_varscan = varscan_file[varscan_file['Position'] == position]['VarAllele'].tolist()
            variants_mutect = mutect_file[mutect_file['POS'] == position]['ALT'].tolist()

            # If there's only one possible variant in both callers, and it's the same in both callers, we add it to Minos
            if len(variants_varscan) == 1 and len(variants_mutect) == 1 and variants_varscan[0] == variants_mutect[0]:
                to_add = [ref_ID[0],str(position), ".", varscan_file[varscan_file['Position'] == position]['Ref'].iloc[0],varscan_file[varscan_file['Position'] == position]['VarAllele'].iloc[0],".",".","VarScan_Mutect2",".","."]
                lines_noheader.append("\t".join(to_add))
            
            # If there's only one possible variant in both callers, but they are different, we add both variants as multiallelic position
            elif len(variants_varscan) == 1 and len(variants_mutect) == 1 and variants_varscan[0] != variants_mutect[0]:
                to_add = [ref_ID[0],str(position), ".", varscan_file[varscan_file['Position'] == position]['Ref'].iloc[0],"{},{}".format(variants_varscan[0], variants_mutect[0]),".",".","VarScan_Mutect2",".","."]
                lines_noheader.append("\t".join(to_add))
            
            # If there's more than one possible variant in any of the callers, we add all the variants as multiallelic position
            elif len(variants_varscan) >= 1 or len(variants_mutect) >= 1:
                all_variants = set(variants_varscan + variants_mutect)
                to_add = [ref_ID[0],str(position), ".", varscan_file[varscan_file['Position'] == position]['Ref'].iloc[0],",".join(all_variants),".",".","VarScan_Mutect2",".","."]
                lines_noheader.append("\t".join(to_add))
    
    # Now we save the Minos output file including the common positions found

    with open("{}.final_sin_wt.vcf_complemented".format(prefix), "w+") as output_file:
        split_data = [line.split('\t') for line in lines_noheader]
        columns = header.lstrip('#').split('\t')
        df = pd.DataFrame(split_data, columns=columns)

        # Now we sort the rows by the position
        df['POS'] = pd.to_numeric(df['POS'])
        df_sorted = df.sort_values(by='POS').reset_index(drop=True)

        # Convert sorted DF to list of strings
        sorted_data_lines = df_sorted.astype(str).apply('\t'.join, axis=1).tolist()
        df_final = lines_preheader + sorted_data_lines

        output_file.write("\n".join(df_final))

    # Annotate complemented Minos VCF with  if H37Rv or MTB_anc reference is used

    cmd = ["java", "-jar", snpEff, "-v", "-hgvs1LetterAa", "-noStats", "-no-downstream", "-no-intergenic", "-no-intron", "-no-upstream", "-noLof",
           "MTB_ancestor", "{}.final_sin_wt.vcf_complemented".format(prefix)]

    if any(item in ['MTB_anc', 'NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome'] for item in ref_ID):
        with open("{}.final_sin_wt_complemented_annotSnpEff.vcf".format(prefix), "w+") as outfh:
            stat = call(cmd, stdout=outfh)
        #UpdateHistory("Anotado con snpEff: java -jar snpEff -v -hgvs1LetterAa -noStats -no-downstream -no-intergenic -no-intron -no-upstream -noLof MTB_ancestor {}.final_sin_wt.vcf_complemented".format(prefix),"custom",prefix)
    
    else:
        sp.run([f"echo {prefix}.final_sin_wt.vcf_complemented was not annotated with SnpEff since the reference is not MTB_anc or H37Rv >> {prefix}.history"], 
        stdout=sp.PIPE, shell=True, universal_newlines=True)

        # Now, we reassign the variable lines_noheader for the next chunk, since up to this point the elements are not sorted 
    with open("{}.final_sin_wt.vcf_complemented".format(prefix), "r+") as filtered_raw_complemented:
        lines_complemented = filtered_raw_complemented.readlines()
        lines_noheader = [x.strip() for x in lines_complemented if "#" not in x]


    with open("{}.minos.raw.tab".format(prefix), "w+") as raw_tab:
        raw_tab.write("#Chrom\tPosition\tRef\tCons\tVarFreq\tCov_allele\tVarAllele\n") # Header

        for line in lines_noheader: # For each position in the Minos output (now complemented with common positions)
            tokens = line.split("\t")
            chrom, position, ref, cons = itemgetter(0,1,3,4)(tokens)

            if "./.:.:0,0:.:0:0,0:0.0:0.0" in line:
                chrom = f"{'-'.join(ref_ID[0:-1])}*"

            varfreq = get_freqs_from_varscan_and_mutect(position, ref_ID) # Get mean frequency from VarScan and Mutect2
            cov_total = get_depths_from_varscan_and_mutect(position, ref_ID) # Get mean depth from VarScan and Mutect2

            if len(cons.split(",")) == 1 and varfreq < 90: # If it's biallelic and the ALT freq is < 90%, we map to IUPAC code
                code_iupac = find_key_with_elements(iupac, [ref, cons])
                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,code_iupac,varfreq,cov_total,cons))
    
            elif len(cons.split(",")) == 1 and varfreq >= 90: # If it's biallelic and the ALT freq is >= 90%, we fix the ALT
                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,cons,varfreq,cov_total,cons))

            elif len(cons.split(",")) > 1 and check_bigger_than_90(varfreq,90): # If it's multiallelic and any ALT freq is >= 90%, we fix that ALT in Cons column but keep all the info in the rest of columns
                variantes = cons.split(",")
                # Get the index of the highest freq value among the possible variants
                max_index = varfreq.index(max(varfreq))
                min_index = varfreq.index(min(varfreq))
                variante_90perc = variantes[max_index] # We get the variant at > 90 % freq
                variante_minoritaria = variantes[min_index] 
                string_varfreq = "{},{}".format(varfreq[max_index], varfreq[min_index])
                string_VarAllele = "{},{}".format(variante_90perc, variante_minoritaria)
                string_depth = "{},{}".format(cov_total[max_index], cov_total[min_index])
                
                # We write the line with the variant > 90% in Cons
                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,variante_90perc,string_varfreq,string_depth,string_VarAllele))
            
            # If there's no variant at > 90% freq, we put ? in Cons and write all the info in parallel order
            elif len(cons.split(",")) > 1 and check_lower_than_90(varfreq,90): 

                variantes = cons.split(",")
                max_index = varfreq.index(max(varfreq))
                min_index = varfreq.index(min(varfreq))
                variante_mayoritaria = variantes[max_index] 
                variante_minoritaria = variantes[min_index] 
                string_varfreq = "{},{}".format(varfreq[max_index], varfreq[min_index])
                string_VarAllele = "{},{}".format(variante_mayoritaria, variante_minoritaria)
                string_depth = "{},{}".format(cov_total[max_index], cov_total[min_index])

                raw_tab.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom,position,ref,"?",string_varfreq,string_depth,string_VarAllele))
                                   
    filter_raw_minos("{}.minos.raw.tab".format(prefix),prefix)
    return 0


def filter_raw_minos(file_name,prefix):
    '''Filters the Minos raw tab output to keep only those positions with depth >= 3 and frequency >= 5%'''
    from operator import itemgetter
    #from .History import UpdateHistory
    # Open the Minos raw tab file
    with open(file_name, "r+") as input_file:
        lines_input = input_file.readlines()

    # Generate list to keep the positions that pass the filter

    posiciones_to_keep = []

    with open("{}.filtered.minos.raw.tab".format(file_name.split(".")[0]), "w+") as output_file:

        headers = [line for line in lines_input if "#" in line]
        lines_noheader = [line for line in lines_input if "#" not in line]
        output_file.write(("").join(headers))

        # Now we filter the positions with depth >= 3 and frequency >= 5%

        for line in lines_noheader: 
            depth, freq = itemgetter(5,4)(line.split("\t"))
            depth = depth.split(",")[0] # If there are multiple depths, we take the highest one, which is the first one
            freq = freq.split(",")[0] # If there are multiple freqs, we take the highest one, which is the first one
            if int(depth) >= 3 and float(freq) >= 5: 
                posiciones_to_keep.append(line)
            else:
                continue

        output_file.write(("").join(posiciones_to_keep)) 
    #UpdateHistory("Minos output filtered keeping variants with depth > 3 and freq > 5", "custom", prefix)
    return 0

  
# EPI and density filtering
def filter_EPI(prefix):
    from .History import UpdateHistory
    '''Keep only those SNPs with depth >= 20 and frequency >= 90%'''
    
    with open("{}.filtered.minos.raw.tab".format(prefix), "r+") as input_file:
        lines = input_file.readlines()
    
    with open("{}.EPI.snp.nodensityfilter".format(prefix), "w+") as output:
        output.write(lines[0]) 
        lines = lines[1:] 

        for line in lines:
            tokens = line.split("\t")
            cons, coverage, var_freq, var_allele = tokens[3], tokens[5], tokens[4], tokens[6]

            if "?" in cons: # If the consensus is unknown, we skip the position
                continue

            elif "," in var_allele and cons in ["A","T","C","G"]: # If it's multiallelic, we check if any allele is fixed and we keep that one
                highest_coverage = coverage.split(",")[0] 
                highest_var_freq = var_freq.split(",")[0] 
                main_variant = tokens[6].split(",")[0]

                if int(highest_coverage) >= 20 and float(highest_var_freq) >= 90:
                    output.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(tokens[0],tokens[1],tokens[2],main_variant,highest_var_freq,
                                                                            highest_coverage,main_variant))
                else:
                    continue

            elif int(coverage) >= 20 and float(var_freq) >= 90:
                output.write(line)
    
    return 0


def densityfilter(prefix, window, dens):
    '''Filter those SNPs that are in high density areas of SNPs'''
    from .History import UpdateHistory
    import pandas as pd
    pos_to_remove = []

    EPI_tab_file = pd.read_csv("{}.EPI.snp.nodensityfilter".format(prefix), sep = '\t')

    # header to inherit
    with open("{}.EPI.snp.nodensityfilter".format(prefix), 'r') as f:
        header = [li for li in f if li.startswith('#')]

    # for each POS
    # how many other POS are less than 'window' pb away from them
    # Is this number greater than dens?
    for p in EPI_tab_file['Position']:
        aux = list(filter(lambda x:
                          x >= (p-window) and x <= (p+window),
                          EPI_tab_file['Position']))
        if len(aux) >= dens:
            pos_to_remove.extend(aux)

    # remove duplicates and sort
    pos_to_remove = list(set(pos_to_remove))
    pos_to_remove.sort()

    # save positions removed by density
    with open("{}.dens_removed_snps.tab".format(prefix), 'w') as rm_file:
        rm_file.write('\t'.join(header))
        EPI_tab_file.loc[EPI_tab_file['Position'].isin(pos_to_remove)].to_csv("{}.dens_removed_snps.tab".format(prefix),mode='a', sep="\t", index=False)


    # save filtered VCF
    final_file = open("{}.EPI.snp.final.annoF".format(prefix), 'w')
        #final_file.write('\t'.join(header))

    EPI_tab_file.loc[~EPI_tab_file['Position'].isin(pos_to_remove)
                 ].to_csv("{}.EPI.snp.final.annoF".format(prefix),
                          mode='a', sep="\t", index=False)
    final_file.close()

    print(len(pos_to_remove), " SNPs removed by density filter: "
          "window = ", window, ", density = ", dens)

    #UpdateHistory("Density filter applied to {}.EPI.snp.vcf file"
    #              " with parameters window = {} and density = "
    #              "{}".format(prefix, window, dens),
    #              "custom", prefix)


def LoadAnnotation(anno_file):
    '''Load the annotation file to store each position to
    be discarded in a list

    This function is meant to load by default
    the H37Rv annotation under
    /data/Databases/MTB_annotation/H37Rv.annotation_new.tsv
    '''
    from .History import UpdateHistory
    discard_positions = []
    with open(anno_file) as infile:
        next(infile)
        for line in infile:
            line = line.strip().split("\t")

            if line[-1] == "DISCARD":
                start, end = int(line[1]), int(line[2])
                for pos in range(start, end + 1):
                    discard_positions.append(str(pos))

    return discard_positions


def correct_line(content):

    iupac = {"M": ["A", "C"],
        "R": ["A", "G"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "Y": ["C", "T"],
        "K": ["G", "T"]}

    amb = iupac[content[3]]

    if (content[2] in amb) and (content[-1] in amb):
        return(content)
    else:
        amb.remove(content[2])
        content[-1] = amb[0]
        return(content)
	
	
def get_DR(prefix):
    '''Perform annotation of DR.snp files and lowcov positions'''
    import os
    from .AnnotatePos import annotateSingle, annotateDouble, annotateTriple

    annotation = "/data/ThePipeline_v3/data/annotation_H37Rv.csv"
    
    ancestor = open("/data/Databases/MTB_ancestor/MTB_ancestor_reference.fasta").read().replace("\n","").replace(">MTB_anc","")
    
    code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"!", "TAG":"!",
    "TGT":"C", "TGC":"C", "TGA":"!", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    genes, names, output, codons = {}, {}, [], []
    
    #first we read the annotation
    with open(annotation) as inAnnot:
        for line in inAnnot:
            line = line.strip().split(",")
            genes[line[0]] = [int(line[-3]),int(line[-2]),line[-1]]
            names[line[0]] = line[1]

    import sys
    #print(genes)

    print("\033[92m\nNow annotating variants from sample {}...\n\033[00m".format(prefix))
    #now we load the DR-snp file and annotate each individual position
    
    with open("{}.filtered.minos.raw.tab".format(prefix)) as infile:
        next(infile)

        for line in infile:
            line = line.strip().split("\t")

            if "," in line[4]: # Modificar esto para que despliegue las líneas con trialélicas en 2 líneas separadas

                freqs = line[4].split(",")
                coverage = line[5].split(",")
                variants = line[6].split(",")
                line_split = [[line[0], line[1], line[2], variants[0], freqs[0], coverage[0], variants[0]],[line[0], line[1], line[2], variants[1], freqs[1], coverage[1], variants[1]]]

                for line_ in line_split:
                # Solo mantenemos las posiciones con una freq mayor de 5% y una profundidad mayor de 6
                    if float(line_[4]) < 5 and float(line_[5]) < 6:
                        continue
                    if line_[3] in ["M", "R", "W", "S", "Y", "K"]:
                        line_ = correct_line(line_)      

                    if len(line_[-1]) == 1:
                        snpid = line_[1]+line_[2]+line_[-1] # POS + REF + ALT
                        change = annotateSingle(snpid, genes, code, ancestor)
                    # print(line)
                        line_.extend(change.strip().split(","))
                        line_.append("triallelic")
                    # print(line)
                        output.append(line_)
            
            else:
                if float(line[4]) < 5 and float(line[5]) < 6:
                    continue
                if line[3] in ["M", "R", "W", "S", "Y", "K"]:
                    line = correct_line(line)      

                if len(line[-1]) == 1:
                    snpid = line[1]+line[2]+line[-1] # POS + REF + ALT
                    change = annotateSingle(snpid, genes, code, ancestor)
                    # print(line)
                    line.extend(change.strip().split(","))
                    # print(line)
                    output.append(line)



    for item in output:
        if "intergenic" not in item:
            codons.append(item[7]+"_"+item[8][1:-1]) # Se guarda la posicion del codon
        else:
            codons.append(item[7])
    
    
    #here we check if there is more than one mutation in the same codon
    double_muts, triple_muts = [], []
    for i in range(0, len(codons)):
        try:
            if codons[i] == codons[i+1]: #we have a multiple mutation

                if codons[i] == codons[i+2]: #this is a triple mutation
                    if "IG" not in codons[i]:
                        triple_muts.append(codons[i])
                        
                else: #it is a double then
                    if "IG" not in codons[i]:
                        double_muts.append(codons[i])

        except IndexError:
            pass
    
    double_muts, triple_muts = list(set(double_muts)), list(set(triple_muts))
    for item in triple_muts:
        if item in double_muts:
            double_muts.remove(item)

    final_output, computed = [], {}
    for i in range(0, len(output)): #now let's correct the preliminary annotation accordingly
        this_codon = output[i][7]+"_"+output[i][8][1:-1]

        if this_codon in double_muts:
            if this_codon not in computed.keys():

                first_mut = output[i][1]+output[i][2]+output[i][6]
                second_mut = output[i+1][1]+output[i+1][2]+output[i+1][6]
                double_change = annotateDouble(first_mut + "+" + second_mut, genes, code, ancestor)

                computed[this_codon] = double_change
            else:
                double_change = computed[this_codon]
                
            new_line = output[i]
            new_line[7:] = double_change.strip().split(",")
            if output[i][1] == output[i+1][1] or output[i][1] == output[i-1][1]:
                new_line.append("TRIALLELIC POSITION")
            else:
                new_line.append("Double mutation")
            final_output.append(new_line)
            
        elif this_codon in triple_muts:
            if this_codon not in computed.keys():
                first_mut = output[i][1]+output[i][2]+output[i][6]
                second_mut = output[i+1][1]+output[i+1][2]+output[i+1][6]
                third_mut = output[i+2][1]+output[i+2][2]+output[i+2][6]

                triple_change = annotateTriple(first_mut + "+" + second_mut + "+" + third_mut, genes, code, ancestor)
                computed[this_codon] = triple_change
            else:
                triple_change = computed[this_codon]
                
            new_line = output[i]
            new_line[7:] = triple_change.strip().split(",")
            if output[i+1][1] == output[i+2][1] or output[i][1] == output[i+1][1] or output[i][1] == output[i-1][1]:
                new_line.append("Double mutation (TRIALLELIC)")
            else:
                new_line.append("Triple mutation")
            final_output.append(new_line)
        else:
            final_output.append(output[i])

    #finally we save the resulting annotated DR.snp.final file
    with open("{}.DR.snp.final".format(prefix), "w") as outfile:
        outfile.write("\t".join(['Chrom', 'Position', 'Ref', 'Cons', 'VarFreq', 'Cov_total', 'VarAllele','Gene','Change','Ann','CodonWT','CodonVAR','Comments\n']))
        for line in final_output:
            outfile.write("\t".join(line)+"\n")
    
    low_cov = {}
    with open("{}.lowcov.tsv".format(prefix)) as infile:
        next(infile)
        for line in infile:
            line = line.strip().split("\t")
            this_gene = ""
            for gene, coords in genes.items():
                if coords[0] <= int(line[0]) <= coords[1]:
                    this_gene = gene
            low_cov[line[0]] = this_gene
    
    with open("{}.lowcov".format(prefix), "w") as outfile:
        outfile.write("Position\tLocus\tGene\n")
        for position, gene in low_cov.items():
            outfile.write("\t".join([position, gene, names[gene]+"\n"]))

def CorrectDR(prefix):
    '''Función para corregir la anotación de las posiciones trialélicas'''

    from .AnnotatePos import annotateSingle, annotateDouble, annotateTriple

    annotation = "/data/ThePipeline_v3/data/annotation_H37Rv.csv"
    
    ancestor = open("/data/Databases/MTB_ancestor/MTB_ancestor_reference.fasta").read().replace("\n","").replace(">MTB_anc","")
    
    code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"!", "TAG":"!",
    "TGT":"C", "TGC":"C", "TGA":"!", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    genes, names, output, codons = {}, {}, [], []
    
    #first we read the annotation
    with open(annotation) as inAnnot:
        for line in inAnnot:
            line = line.strip().split(",")
            genes[line[0]] = [int(line[-3]),int(line[-2]),line[-1]]
            names[line[0]] = line[1]

    with open("{}.DR.snp.final".format(prefix), "r+") as DR_file:
        lines = DR_file.readlines()

    with open("{}.DR.snp.final".format(prefix), "w+") as output_DR_file:
        output_DR_file.write(lines[0])

        to_skip = []

        for i in range(1, len(lines)):
            line = lines[i].strip().split("\t")

            if line[1] in to_skip:
                continue

            # If the triallelic position is not part of a double mutation with another position 
            # in the same codon, we re‑annotate it independently for each allele

            if "TRIALLELIC POSITION" in line: 
                new_line = line[0:7]
                snpid = new_line[1]+new_line[2]+new_line[-1] # POS + REF + ALT
                change = annotateSingle(snpid, genes, code, ancestor)
                new_line.extend(change.strip().split(","))
                new_line.append("TRIALLELIC POSITION")
                output_DR_file.write("\t".join(new_line) + "\n")
            
            elif "Double mutation (TRIALLELIC)" in line:

                first_mut = line[1]+line[2]+line[6]
                second_mut = lines[i+1].split("\t")[1]+lines[i+1].split("\t")[2]+lines[i+1].split("\t")[6]
                first_double_change = annotateDouble(first_mut + "+" + second_mut, genes, code, ancestor)
                second_mut = lines[i+2].split("\t")[1]+lines[i+2].split("\t")[2]+lines[i+2].split("\t")[6]
                second_double_change = annotateDouble(first_mut + "+" + second_mut, genes, code, ancestor)

                new_line_first = line[0:7]
                new_line_first.extend(first_double_change.strip().split(","))
                new_line_first.append("Double mutation (TRIALLELIC)")
                new_line_second = line[0:7]
                new_line_second.extend(second_double_change.strip().split(","))
                new_line_second.append("Double mutation (TRIALLELIC)")

                new_line_third = lines[i+1].split("\t")[0:7]
                new_line_third.extend(first_double_change.strip().split(","))
                new_line_third.append("Double mutation (TRIALLELIC)")

                new_line_fourth = lines[i+2].split("\t")[0:7]
                new_line_fourth.extend(second_double_change.strip().split(","))
                new_line_fourth.append("Double mutation (TRIALLELIC)")

                to_skip.append(line[1])
                to_skip.append(lines[i+1].split("\t")[1])

                output_DR_file.write("\t".join(new_line_first) + "\n" + "\t".join(new_line_third) + "\n" + "\t".join(new_line_second) + "\n" + "\t".join(new_line_fourth) + "\n")
            
            else:
                output_DR_file.write("\t".join(line)+"\n")

    
    
def Calling(args):
    from .Repository import Programs, Data
    import sys
    import os
    import subprocess as sp
    import shutil
    # Check package dependencies
    try:
        import vcf
    except ModuleNotFoundError:
        sys.exit("\033[91mPIPELINE ERROR: The VCF module (PyVCF) is not installed "
                 "and is necessary for the 'Calling' process\033[0m")
    try:
        import pandas
    except ModuleNotFoundError:
        sys.exit("\033[91mPIPELINE ERROR: The PANDAS module is not installed "
                 "and is necessary for the 'Calling' process\033[0m")

    programs = Programs()
    data = Data()

    varscan = programs["varscan"]
    minos = programs["minos"]
    gatk = programs["gatk"]
    snpEff = programs["snpEff"]
    samtools = programs["samtools"]
    reference = args.reference
    genomeCoverageBed = programs["genomeCoverageBed"]

    # Check that mapping files end with .bam/.cram
    check_extension(args.ext)

    if not reference:
        reference = data["reference"]
    
    ref_ID = get_reference_id(args.prefix, args.ext, samtools)

    # SNP calling with VarScan2
    print("\033[92m\nRunning VarScan on {}\n\033[00m".format(args.prefix))
    VarScan(reference, args.prefix, varscan, samtools, args.ext, ref_ID)

    # SNP calling with Mutect2
    print("\033[92m\nRunning Mutect2 on {}\n\033[00m".format(args.prefix))
    Mutect2(reference, args.prefix, gatk, samtools, genomeCoverageBed,
            args.ext, args.threads, args.min_depth, args.min_qual,
            args.filter_freq)

    print("\033[92m\nAnnotation-filtering VarScan and Mutect2 outputs of sample {}\n\033[00m".format(args.prefix))
    annoF_varscan_mutect(args.prefix)

    # Multiallelic filter from Mutect2 snp output
    filter_multiallelic_from_mutect2_snp(args.prefix)

    #let's call the wild type positions from the gVCF file
    print("\033[92m\nCalling WT for {}\n\033[00m".format(args.prefix))
    callWT(args.prefix)

    # Minos calling

    #we now copy the reference sequence for Minos analysis to avoid
    #unnecessary singularity path bindings, and run the program inside its container
    os.popen('cp ' + reference + ' .')
    print("\033[92m\nPerforming variant call adjudication with Minos for sample {}...\n\033[00m".format(args.prefix))
    Minos(reference,minos,snpEff,args.prefix, args.single_end, ref_ID)
    minos_raw_vcf_to_tab(args.prefix, ref_ID, snpEff)

    # EPI filtering
    print("\033[92m\nObtaining {}.EPI.snp.final.annoF\n\033[00m".format(args.prefix))
    filter_EPI(args.prefix)
        
    # Apply density filter
    if not args.skip_dens_filt:  # by density
        print("\033[92m\nFiltering {}.EPI.snp.final.annoF by density\n\033[00m".format(args.prefix))
        densityfilter(args.prefix, args.window, args.density)


    # DR extraction and correction
    print("\033[92m\nObtaining {}.DR.snp.final\n\033[00m".format(args.prefix, args.prefix))
    get_DR(args.prefix)
    CorrectDR(args.prefix)
    
    print("\033[92m\nRemoving and renaming files\n\033[00m")
    #a final cleanup and renaming of intermediate and/or useless files
    os.remove("{}.EPI.snp.nodensityfilter".format(args.prefix))
    os.rename("{}.filtered.minos.raw.tab".format(args.prefix), "{}.snp.minos".format(args.prefix))
    os.remove("{}.final_sin_wt.vcf".format(args.prefix))
    os.remove("{}.lowcov.tsv".format(args.prefix))
    os.remove("{}.minos.raw.tab".format(args.prefix))
    os.remove("{}.parsed.vcf".format(args.prefix))
    os.rename("{}.remade.snp.vcf".format(args.prefix),"{}.snp.mutect".format(args.prefix))
    os.rename("{}.snp".format(args.prefix),"{}.snp.varscan".format(args.prefix))
    os.remove("{}.snp.vcf".format(args.prefix))
    os.remove("{}.vcf".format(args.prefix))
    os.remove("{}_mapq60.sort.bam".format(args.prefix))
    os.remove("{}_no_orientation.vcf".format(args.prefix))
    os.remove("{}_unfiltered.vcf".format(args.prefix))
    os.remove("{}_unfiltered.vcf.stats".format(args.prefix))
    os.remove("{}.final_sin_wt.vcf_complemented".format(args.prefix))
    os.remove("{}.f1r2.tar.gz".format(args.prefix))
    os.remove("{}.vcf.filteringStats.tsv".format(args.prefix))
    os.remove("{}.vcf.idx".format(args.prefix))
    os.remove("read-orientation-model_{}.tar.gz".format(args.prefix))
    # If not -kna, remove original SNP files from VarScan and Mutect2 and just keep files filtered by annotation
    if not args.keep_not_annof:
        os.remove("{}.parsed.vcf.original_no_annoF".format(args.prefix))
        os.remove("{}.snp.original_no_annoF".format(args.prefix))
        os.remove("{}.snp.vcf.original_no_annoF".format(args.prefix))


    # Mutect2 output VCF file is converted to tab format by default if -kmvcf parameter is not used. If it is, we keep both files
    mutect2_vcf_to_tab(args.prefix)
    if not args.kmvcf:
        os.remove("{}.snp.mutect".format(args.prefix))


    # If the reference used is not based in H37Rv (i.e. MTB_anc or H37Rv), gene annotation won't be correct. Thus, we remove the annotation from the annotated files.
    if not any(item in ['MTB_anc', 'NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome'] for item in ref_ID):
        sp.run([f"cut -f1-7 {args.prefix}.DR.snp.final > {args.prefix}.DR.temp"], stdout=sp.PIPE, shell=True, universal_newlines=True)
        sp.run([f"mv {args.prefix}.DR.temp {args.prefix}.DR.snp.final"], stdout=sp.PIPE, shell=True, universal_newlines=True)
        sp.run([f"cut -f 1 {args.prefix}.lowcov > {args.prefix}.lowcov.temp"], stdout=sp.PIPE, shell=True, universal_newlines=True)
        sp.run([f"mv {args.prefix}.lowcov.temp {args.prefix}.lowcov"], stdout=sp.PIPE, shell=True, universal_newlines=True)
        sp.run([f"echo {args.prefix}.DR.snp.final and {args.prefix}.lowcov annotations were removed since the used reference was not based on H37Rv >> {args.prefix}.history"],
        stdout=sp.PIPE, shell=True, universal_newlines=True)


    if not args.keep_gvcf:
        os.remove("{}.gvcf".format(args.prefix))
        os.remove("{}.gvcf.stats".format(args.prefix))
    else:
        print("\033[92m\nCompressing {}.gvcf\n\033[00m".format(args.prefix))
        sp.run("gzip {}.gvcf".format(args.prefix), shell=True)

    if not args.keep_bam:
        # If .cram file not present, convert bam into cram and then remove bam
        if "{}.cram".format(args.prefix) not in os.listdir():
            cmd_toCRAM = [samtools, "view", "--threads", args.threads,
            "-T", reference,
            "-C", "-o", "{}.cram".format(args.prefix),
            "{}{}".format(args.prefix, args.ext)]
            sp.call(cmd_toCRAM)
            os.remove("{}.sort.bam".format(args.prefix))
            os.remove("{}.sort.bam.bai".format(args.prefix))
        else:
            os.remove("{}.sort.bam".format(args.prefix))
            os.remove("{}.sort.bam.bai".format(args.prefix))
    
    return 0
