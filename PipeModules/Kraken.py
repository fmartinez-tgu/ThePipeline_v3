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

from sys import stdout


def KrakenClassify(fastq, krakenDB, prefix, paired, threads, compressed,
                   kraken, preload):
    '''Run kraken to classify fastq reads using krakenDB and store results
    in prefix.kraken file'''
    from subprocess import call
    import sys
    from .History import UpdateHistory

    cmd = [kraken, "--db", krakenDB,
           "--fastq-input", "--threads", threads]
    if compressed:
        cmd.append("--gzip-compressed")
    if preload:
        cmd.append("--preload")
    if paired:
        cmd.append("--paired")

    for fq in fastq:
        cmd.append(fq)

    with open("{}.kraken".format(prefix), "w") as outfh:
        stat = call(cmd, stdout=outfh)

    if stat != 0:
        sys.exit("\033[91mERROR: Pipeline stopped when"
                 " classifying reads with Kraken!\n\033[0m")

    UpdateHistory(cmd, "kraken", prefix)
    return 0


def KrakenTranslate(krakenDB, prefix, kraken_translate):
    ''' Translate .kraken file to .labels file '''
    from subprocess import call
    import sys
    from .History import UpdateHistory

    cmd = [kraken_translate, "{}.kraken".format(prefix), "--db", krakenDB]
    with open("{}.labels".format(prefix), "w") as outfh:
        stat = call(cmd, stdout=outfh)

    if stat != 0:
        sys.exit("\033[91mERROR: Pipeline stopped "
                 "when translating kraken file!\n\033[0m")

    UpdateHistory(cmd, "kraken", prefix)
    return 0


def MakeReadList(mstring, prefix):
    '''Call grep to select those reads that match the string
    This string is Mycobacterium tuberculosis by default'''
    import subprocess
    import sys
    from .History import UpdateHistory

    cmd = ["/usr/bin/fgrep", mstring, "{}.labels".format(prefix), "|",
           "cut", "-f1", ">", "{}.filtered.readlist".format(prefix)]
    ps = subprocess.Popen(("/usr/bin/fgrep", mstring,
                           "{}.labels".format(prefix)),
                          stdout=subprocess.PIPE)
    output = subprocess.check_output(("cut", "-f1"), stdin=ps.stdout)
    ps.wait()

    with open("{}.filtered.readlist".format(prefix), "w") as outfh:
        # IMPORTANT! While typing 'output' prints data in pyhton2.7, in
        # python3+ you must 'decode' the output of the
        # subprocess. Otherwise, python3 interprets it as 'byte'
        # instead of 'str'
        # check https://docs.python.org/3.0/whatsnew/3.0.html#text-vs-data-instead-of-unicode-vs-8-bit
        outfh.write(output.decode(sys.stdout.encoding))

    UpdateHistory(cmd, "kraken", prefix)
    return 0


def PickReads(fastq, prefix, seqtk, pigz, threads, strand):
    '''Write reads in readlist to output'''
    from subprocess import call
    import sys
    from .History import UpdateHistory

    readlist = "{}.filtered.readlist".format(prefix)
    cmd = [seqtk, "subseq", fastq, readlist, "|",
           "pigz -p ", threads, ">",
           "{}{}.filtered.fastq".format(prefix, strand)]
    pseqtk = [seqtk, "subseq", fastq, readlist]
    with open("{}{}.filtered.fastq".format(prefix, strand), "w") as outfh:
        stat = call(pseqtk, stdout=outfh)
        if stat != 0:
            sys.exit("\033[91mERROR: Pipeline stopped when"
                     " selecting read from readlist! "
                     "Are you trying to run the --filter option before "
                     "running the --classify option?\n\033[0m")

    # parallel gzip of filtered fastq
    pigz_ps = [pigz, "-p", threads, "-f",
               "{}{}.filtered.fastq".format(prefix, strand)]
    stat = call(pigz_ps, stdout=stdout)

    UpdateHistory(cmd, "seqtk", prefix)
    return 0


def clean_files(prefix):
    '''Clean .kraken, .labels and .readlist files once .filtered.readlist file
    has been created'''
    import os

    os.remove("{}.labels".format(prefix))
    os.remove("{}.filtered.readlist".format(prefix))
    os.remove("{}.kraken".format(prefix))

    return 0


def MakeReport(prefix, krakenDB, version=True):
    '''This function calls kraken-report if specified'''
    import subprocess as sp
    import sys
    from .History import UpdateHistory

    cmd = ["/data/ThePipeline_programs/Kraken/kraken-report",
           "--db", krakenDB, "{}.kraken".format(prefix)]
    with open("{}.kraken.report".format(prefix), "w") as outfh:
        stat = sp.call(cmd, stdout=outfh)

    if stat != 0:
        sys.exit("\033[91mERROR: Pipeline stopped"
                 " when making the kraken report!\n\033[0m")

    UpdateHistory(cmd, "kraken", prefix)
    return 0


def AssessContam(prefix):

    '''Write genus and species tables'''
    genus_table, species_table = ReportParser(prefix)

    with open("{}.genus.contaminants".format(prefix), "w") as outfile:
        outfile.write("Sample\tClassification\tReads\tPercentage\n")
        for row in genus_table:
            sample, taxonomy, reads, perc = row
            outfile.write("{}\t{}\t{}\t{}\n".format(sample, taxonomy, reads,
                          perc))

    with open("{}.species.contaminants".format(prefix), "w") as outfile:
        outfile.write("Sample\tClassification\tReads\tPercentage\n")
        for row in species_table:
            sample, taxonomy, reads, perc = row
            outfile.write("{}\t{}\t{}\t{}\n".format(sample, taxonomy, reads,
                          perc))


def ReportParser(prefix):
    '''Parse a kraken-report and produce two tables (genus and species) with
    contaminants and its frequency'''
    MTB_species = ["Mycobacterium tuberculosis",
                   "Mycobacterium africanum",
                   "Mycobacterium bovis",
                   "Mycobacterium canettii",
                   "Mycobacterium microti",
                   "Mycobacterium caprae",
                   "Mycobacterium pinnipedii"]
    # Tables are stored in a list so it can be
    # The ordered according to the frequency of
    # contaminants
    genus_table = []
    species_table = []

    # We cannot assume there is going to be mycobacterium
    # so initialize one variable in case myc_reads and myc_freq
    # is never assigned
    myc_reads = 0
    myc_freq = 0.0
    mtb_reads = 0
    mtb_freq = 0

    with open("{}.kraken.report".format(prefix)) as infile:
        unclass_line = infile.readline()
        unclass_line = unclass_line.rstrip()
        freq, rcov, rassign, level, taxid, taxonomy = unclass_line.split("\t")
        row = [prefix, "Unclassified", int(rcov), float(freq)]
        genus_table.append(row)
        species_table.append(row)

        root_line = infile.readline()
        root_line = root_line.rstrip()
        freq, rcov, rassign, level, taxid, taxonomy = root_line.split("\t")

        for line in infile:
            line = line.rstrip()
            freq, rcov, rassign, level, taxid, taxonomy = line.split("\t")
            taxonomy = taxonomy.strip()  # Remove white spaces
            freq = float(freq)
            rcov = int(rcov)
            if level == "G":
                if taxonomy == "Mycobacterium":
                    myc_reads = rcov
                    myc_freq = freq
                else:
                    row = [prefix, taxonomy, rcov, freq]
                    genus_table.append(row)
            elif level == "S":
                # We're gonna treat MTB and nonMTB species differently
                # as all MTB species are gonna be classified as MTBC
                if taxonomy not in MTB_species:
                    row = [prefix, taxonomy, rcov, freq]
                    species_table.append(row)
            elif taxonomy == "Mycobacterium tuberculosis complex":
                mtb_reads = rcov
                mtb_freq = freq
                row = [prefix, taxonomy, rcov, freq]
                species_table.append(row)
                genus_table.append(row)

        #  Take away MTB reads from Mycobacterium genus

        myc_reads -= mtb_reads
        myc_freq -= mtb_freq
        myc_row = [prefix, "Non-TB-Mycobacterium", myc_reads, myc_freq]
        genus_table.append(myc_row)

    # Sort tables by decreasing percentaje
    genus_table = sorted(genus_table, key=lambda x: x[3], reverse=True)
    species_table = sorted(species_table, key=lambda x: x[3], reverse=True)

    return genus_table, species_table

def CountReads(prefix):
    '''Count total reads analyzed and MTBC-classified reads'''
    from subprocess import check_output
    import sys

    try:
        with open("{}.nfilter".format(prefix), "w") as outfile:
            analyzed = check_output(["wc", "-l", "{}.kraken".format(prefix)])
            analyzed = analyzed.split()[0]
            MTBC = check_output(["wc", "-l", "{}.filtered.readlist".format(prefix)])
            MTBC = MTBC.split()[0]

            outfile.write("TOTAL:{}\n".format(analyzed.decode(sys.stdout.encoding)))#cambio Mariana
            outfile.write("MTBC:{}\n".format(MTBC.decode(sys.stdout.encoding)))#cambio Mariana
            return 0

    except:
        sys.exit("Either .kraken or .filtered.readlist files do not exist")

    assert False



def Kraken(args):
    ''' Launch all function to classify with kraken, translate labels and
    pick reads matching mstring'''
    from .Repository import Programs, Data
    import sys

    programs = Programs()
    data = Data()

    if not args.krakenDB:
        args.krakenDB = data["krakenDB"]

    kraken = programs["kraken"]
    kraken_translate = programs["kraken_translate"]
    seqtk = programs["seqtk"]
    pigz = programs["pigz"]

    if not args.classify and not args.filterf and not args.report:
        if not args.report:
            sys.exit("\033[91m\nERROR: ThePipeline kraken must be"
                     " called with at least one"
                     " of the following options:\n --classify/--filter/"
                     " --report. All can be used at the same time.\n\033[0m")
    if args.classify:
        if not args.fastq:
            sys.exit("\033[91m\nERROR: ThePipeline3 kraken"
                     " with the --classify option needs"
                     " input fastq files! "
                     "\nPass them using the -f argument.\n\033[0m")
        KrakenClassify(args.fastq, args.krakenDB, args.prefix, args.paired,
                       args.threads, args.compressed, kraken, args.preload)
        KrakenTranslate(args.krakenDB, args.prefix, kraken_translate)
        MakeReadList(args.matching, args.prefix)
        CountReads(args.prefix)
        MakeReport(args.prefix, args.krakenDB)

    if args.filterf:
        if not args.fastq:
            sys.exit("\033[91m\nThePipeline3 kraken with"
                     " the --filter option needs"
                     " input fastq files! "
                     "\nPass them using the -f argument.\n\033[0m")
        elif len(args.fastq) == 1:
            PickReads(args.fastq[0], args.prefix, seqtk, pigz, args.threads,
                      strand="")
        elif len(args.fastq) == 2:
            PickReads(args.fastq[0], args.prefix, seqtk, pigz, args.threads,
                      strand=".P1")
            PickReads(args.fastq[1], args.prefix, seqtk, pigz, args.threads,
                      strand=".P2")
        else:
            assert False  # change this and manage exception

    if args.report:
        AssessContam(args.prefix)

    if args.noclean:
        return 0
    else:
        try:
            clean_files(args.prefix)
            return 0
        except OSError:
            # We're trying to remove kraken and labels files
            # and moving clean.fastq files
            # but maybe they don't exist
            return 0
