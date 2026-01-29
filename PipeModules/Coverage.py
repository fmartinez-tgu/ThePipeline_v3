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

from zipfile import Path


def CoverageBAM(prefix, reference):
    ''' Calculate coverage '''
    from subprocess import call
    from .Repository import Programs
    import sys

    # genomeCoverageBed is a script that is under bedtools/bin folder
    programs = Programs()
    genomeCoverageBed = programs["genomeCoverageBed"]

    bamfile = "{}.sort.bam".format(prefix)
    covfile = "{}.coverage".format(prefix)
    with open(covfile, "w") as outfh:
        stat = call([genomeCoverageBed, "-ibam", bamfile, "-d", "-g",
                    reference], stdout=outfh)
    if stat != 0:
        sys.exit("\033[91mERROR: Pipeline stopped when"
                 " calculating each bam coverage!\n\033[0m")

    return 0


def CoverageCRAM(prefix, reference):
    ''' Calculate coverage '''
    import subprocess as sp
    from .Repository import Programs

    # genomeCoverageBed is a script that is under bedtools/bin folder
    programs = Programs()
    genomeCoverageBed = programs["genomeCoverageBed"]
    samtools = programs["samtools"]

    cramfile = "{}.cram".format(prefix)
    covfile = "{}.coverage".format(prefix)
    with open(covfile, "w") as outfh:
        view = sp.Popen([samtools, "view", "-b", cramfile], stdout=sp.PIPE)

        genomeCov = sp.Popen([genomeCoverageBed, "-ibam", "stdin", "-d",
                             "-g", reference], stdin=view.stdout,
                             stdout=outfh)
        view.stdout.close()
        out, err = genomeCov.communicate()

    return 0


def MeanCoverage(prefix, depth4cov):
    ''' Calculate the mean coverage per pb in a sample from
    prefix.coverage file '''

    positions = 0.0
    coverage = 0
    covered_positions = 0.0
    # Add cov list to sort and calc median
    covlist = []
    with open("{}.coverage".format(prefix)) as infile:
        for line in infile:
            if line == "":
                positions = 0.0
                break
            ref, pos, cov = line.rstrip().split()
            cov = int(cov)
            positions += 1.0
            coverage += cov
            if cov >= depth4cov:
                covered_positions += 1.0
            covlist.append(int(cov))

    if positions > 0:
        # Calc mean
        mean_cov = coverage / positions
        # Calc median
        # Recall that position X in a list is accesed with X-1
        # so for example the 6th element is accessed as covlist[5]
        covlist.sort()
        if len(covlist) % 2 == 0:
            m1 = int(len(covlist) / 2)
            cov1 = covlist[m1 - 1]
            cov2 = covlist[m1]
            median_cov = (cov1 + cov2) / 2.0
        else:
            m1 = len(covlist) // 2
            # not neccessary to sum 1 as the m1 +1 element is
            # accessed in fact with just m1
            median_cov = covlist[m1]

        genome_coverage = covered_positions / positions
        return (mean_cov, median_cov, genome_coverage)
    else:
        return (0, 0, 0)


def MeanSampleCoverage(prefix, keepcoverage, depth4cov):
    import os

    meancov, mediancov, genome_coverage = MeanCoverage(prefix, depth4cov)
    with open("{}.meancov".format(prefix), "w") as outfh:
        outfh.write("{}_mean_depth\t{}\n".format(prefix, meancov))
        outfh.write("{}_median_depth\t{}\n".format(prefix, mediancov))
        outfh.write("{}_genome_coverage\t{}\n".format(prefix, genome_coverage))
        outfh.write("\nNote: Minimum depth for considering a base covered "
                    "was {}\n".format(depth4cov))

    # Remove .coverage files
    if not keepcoverage:
        os.remove("{}.coverage".format(prefix))

    return meancov, mediancov, genome_coverage


def FiltByCov(prefix, mean, median, cov, minmean, minmedian, mincov):
    '''If a sample do not pass coverage/depth thresholds, move all its files
    to NoPassCov folder'''
    import os
    import shutil
    from glob import glob

    if mean < minmean or median < minmedian or cov < mincov:
        # Create NoPassCov if it does not exist
        try:
            os.mkdir("zNoPassCov")
        except OSError:
            pass

        files = glob("{}*".format(prefix))
        for file in files:
            shutil.move(file, "./zNoPassCov/{}".format(file))

    return 0


def CalcCoverage(args):

    from .Repository import Data
    import os

    data = Data()

    if not args.reference:
        args.reference = data["reference"]

    # Let's check if prefix.coverage is already present, since it's called in the Calling module

    file_path_coverage = f"{args.prefix}.coverage"

    if not os.path.isfile(file_path_coverage):
        if args.extension == "bam":
            # Calc coverage given a sort.bam or sort.cram
            CoverageBAM(args.prefix, args.reference)
        elif args.extension == "cram":
            CoverageCRAM(args.prefix, args.reference)
        else:
            assert False, "BAM or CRAM extension must be specified."

    # Create the .meancov file and get its values
    mean, median, cov = MeanSampleCoverage(args.prefix, args.keepcoverage,
                                           args.depth4cov)

    # Filter by Coverage in case it is specified
    if args.filter:
        FiltByCov(args.prefix, mean, median, cov,
                  args.minmean, args.minmedian, args.mincov)

    return 0
