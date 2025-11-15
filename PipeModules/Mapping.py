#! /usr/bin/env python3.7
# Copyrigth (C) 2025 Alvaro Chiner Oms & Francisco Jose Martinez

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


def bwa_map(fastq, reference, threads, index, prefix, mapq_cutoff, cram_compress):
    ''' Use bwa to map input fastq/s to reference. It is always expected that
    all programs this pipeline depends on are stored under Programs/

    When no reference is provided, the MTBC ancestor reference (path defined in
    data_path) is used by default.

    Reference is indexed thorugh bwa index if -i/--index option is used or
    directly used as reference otherwise

    Exit status of each command is stored in stat and pipeline stops if one
    of them fails

    stat == 0 means a command succeded
    stat != 0 means a command failed

    '''
    import subprocess as sp
    import sys
    from .Repository import Programs
    from .History import UpdateHistory

    programs = Programs()
    bwa = programs["bwa"]
    samtools = programs["samtools"]

    if index:  # if -i/--index option provided index reference
        stat = sp.call([bwa, "index", "-a", "bwtsw", reference])
        # Check that command succeeded
        if stat != 0:
            sys.exit("\033[91mPIPELINE ERROR at sample {}:"
                     " Pipeline stopped at"
                     " bwa index point!\n\033[0m".format(prefix))

    # Read Group to pass to bwa to make it compatible with GATK
    rgroup = '@RG\\tID:' + prefix + '\\tSM:' + prefix
    rgroup += '_sm\\tPU:' + prefix + '_pu\\tLB:' + prefix + '_lb'

    cmd_bwa = [bwa, "mem", "-t", threads, "-R", rgroup, reference]

    # Add an awk line to filter alignments with MAPQ < mapq_cutoff

    cmd_awk = ['awk', ''.join(['$1 ~ /^@/ || $5 >= ', str(mapq_cutoff)])]
    cmd_sam = [samtools, "view", "-bt", reference, "-", "--threads", threads]
    cmd_sort = [samtools, "sort",  "-o", "{}.sort.bam".format(prefix),
                "--threads", threads]

    for fq in fastq:
        cmd_bwa.append(fq)
    mem = sp.Popen(cmd_bwa, stdout=sp.PIPE)
    if mapq_cutoff > 0:
        awk_MAPQ = sp.Popen(cmd_awk, stdin=mem.stdout, stdout=sp.PIPE)
        mem.stdout.close()
        view = sp.Popen(cmd_sam, stdin=awk_MAPQ.stdout, stdout=sp.PIPE)
        awk_MAPQ.stdout.close()
        sort = sp.Popen(cmd_sort, stdin=view.stdout)
        view.stdout.close()
        out, err = sort.communicate()  # What are out and err used for?
    else:
        view = sp.Popen(cmd_sam, stdin=mem.stdout, stdout=sp.PIPE)
        mem.stdout.close()
        sort = sp.Popen(cmd_sort, stdin=view.stdout)
        view.stdout.close()
        out, err = sort.communicate()  # What are out and err used for?

    UpdateHistory(cmd_bwa, "bwa", prefix)
    if mapq_cutoff > 0:
        UpdateHistory(cmd_awk, "awk", prefix)
    UpdateHistory(cmd_sam, "samtools", prefix)
    UpdateHistory(cmd_sort, "samtools", prefix)

    return 0


def RemoveDuplicates(prefix, keepduplicatesbam):
    '''This function builds a CMD to mark and remove duplicates from .sort.bam
    files using picard. Default parameters are:

    ASSUME_SORTED=true
    REMOVE_DUPLICATES=true
    VALIDATION_STRINGENCY=LENIENT
    '''
    import subprocess as sp
    import sys
    import os
    from .Repository import Programs
    from .History import UpdateHistory

    programs = Programs()
    picard = programs["picard"]

    CMD = ["java", "-Xmx8g", "-jar", picard, "MarkDuplicates"]
    CMD.append("I={}.sort.bam".format(prefix))
    CMD.append("O={}.nd.sort.bam".format(prefix))
    CMD.append("M={}.dup.metrix".format(prefix))
    CMD.append("ASSUME_SORTED=true")
    CMD.append("REMOVE_DUPLICATES=true")
    CMD.append("VALIDATION_STRINGENCY=LENIENT")

    stat = sp.call(CMD)
    if stat != 0:
        sys.exit("\033[91mPIPELINE ERROR at sample {}: Pipeline stopped at"
                 " picard when removing duplicates!\n\033[0m".format(prefix))
    UpdateHistory(CMD, "picard", prefix)

    # Now mv nd.sort.bam to .sort.bam extension so it is compatible with
    # subsequent analysis steps in the pipeline
    if not keepduplicatesbam:
        os.rename("{}.nd.sort.bam".format(prefix),
                  "{}.sort.bam".format(prefix))
    else:
        os.rename("{}.sort.bam".format(prefix),
                  "{}.dup.sort.bam".format(prefix))
        os.rename("{}.nd.sort.bam".format(prefix),
                  "{}.sort.bam".format(prefix))

    return 0

def hardclipping(prefix, reference, samtools, samclip_h):
    ''' This function performs hardclipping filtering of reads using samclip_h tool'''

    import subprocess as sp
    import os
    cmd = f"{samtools} view -h {prefix}.sort.bam | {samclip_h} --ref {reference} | {samtools} sort -o {prefix}_hardclipping.sort.bam"
    sp.run(cmd, shell=True, check=True)
    os.remove(f"{prefix}.sort.bam")
    os.rename(f"{prefix}_hardclipping.sort.bam",f"{prefix}.sort.bam")


def Mapping(args):
    '''Call previous functions to perform a complete subpipeline'''
    import sys
    import os
    import subprocess as sp
    from .Repository import Programs, Data
    from .History import UpdateHistory

    data = Data()
    programs = Programs()
    samtools = programs["samtools"]
    qualimap = programs["qualimap"]
    samclip_h = programs["samclip_h"]

    # Asign reference
    if not args.reference:
        args.reference = data["reference"]

    # Check that only one or two fastqs have been provided
    if len(args.fastq) != 2 and len(args.fastq) != 1:
        sys.exit("\033[91mPIPELINE ERROR at sample {}:"
                 " Up to two fastq "
                 "files can be provided for mapping.\n\033[0m".format(args.prefix))
    if args.cram_compress:
        cram_compress = True
    else:
        cram_compress = False

    bwa_map(args.fastq, args.reference, args.threads, args.index, args.prefix,
            args.mapq_cutoff, args.cram_compress)

    if args.nodedup is False:
        RemoveDuplicates(args.prefix, args.keepduplicatesbam)
    
    if not args.nhc:
        hardclipping(args.prefix, args.reference, samtools, samclip_h)


    # Execute QualiMap
    # You need to unset DISPLAY to avoid errors
    try:
        sp.run("unset DISPLAY", shell=True,
               capture_output=True)
        cmd_qualimap = [qualimap, "bamqc", "-bam",
                        "{}.sort.bam".format(args.prefix)]
        # reset DISPLAY so you can open graphic envinronments
        sp.run("export $DISPLAY='{}'".format(os.environ['DISPLAY']),
               shell=True, capture_output=True)

        stat = sp.call(cmd_qualimap)
        if stat != 0:
            sys.stdout.write("\033[93mPIPELINE WARNING at sample {}: Error "
                             "at QualiMap"
                             " execution!\n\033[0m".format(args.prefix))
    except:
        sys.stdout.write("\033[93mPIPELINE WARNING at sample {}: Error"
                         " unsetting DISPLAY. Have you set the -X option"
                         " in the ssh connection? QualiMap has not"
                         " been executed. \n\033[0m".format(args.prefix))

    # Optionally, transform bam to cram
    if cram_compress:
        cmd_toCram = [samtools, "view", "-T", args.reference, "-C", "-o",
                  "{}.cram".format(args.prefix),
                  "{}.sort.bam".format(args.prefix)]

        stat = sp.call(cmd_toCram)
        if stat != 0:
            sys.exit("\033[91mPIPELINE ERROR at sample {}:"
                 " Pipeline stopped at samtools "
                 "when transforming BAM to CRAM!\n\033[0m".format(args.prefix))
        else:
            os.remove("{}.sort.bam".format(args.prefix))

        UpdateHistory(cmd_toCram, "samtools", args.prefix)

    return 0
