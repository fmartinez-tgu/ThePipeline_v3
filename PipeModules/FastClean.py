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

def ReadConfig(cfgfile):
    '''This function read fastP paramters from a config file. If no config file
    is provided, then it reads from ThePipeline_v3/data/Configs/
    It returns a list of parameters ready
    to build the command line (CMD) to call fastP. In the default config file
    these are the specified parameters:
    --cut_by_quality3
    --cut_window_size=10
    --cut_mean_quality=20
    --length_required=50
    --correction
    '''

    import os
    import sys

    parameters = []

    if cfgfile == "default_config":
        config = "{}/data/Configs/fastp_default.config".format(
            os.path.split(
                os.path.dirname(os.path.abspath(__file__)))[0])
        infile = open(config)
    else:
        try:
            infile = open(cfgfile)
        except IOError:
            sys.exit("\033[91mERROR: Configuration file not found\033[0m")

    for line in infile:
        if not line.startswith("#"):
            parameter = line.rstrip()
            parameters.append(parameter)

    infile.close()

    return parameters


def BuildCMD(args):
    ''' This function builds a list with the command line (CMD) that must be
    passed to "call" in order to call fastp'''
    import sys
    from .Repository import Programs

    # Get fastp PATH.
    programs = Programs()
    fastp = programs["fastp"]

    CMD = [fastp]

    # Add input and output files to CMD
    if len(args.fastq) == 1:
        CMD.append("--in1={}".format(args.fastq[0]))
        CMD.append("--out1={}.clean.fastq.gz".format(args.prefix))
    elif len(args.fastq) == 2:
        CMD.append("--in1={}".format(args.fastq[0]))
        CMD.append("--out1={}.P1.clean.fastq.gz".format(args.prefix))
        CMD.append("--in2={}".format(args.fastq[1]))
        CMD.append("--out2={}.P2.clean.fastq.gz".format(args.prefix))
    else:
        sys.exit("\033[91mERROR: You must provide one or two fastq files."
                 " Not less; not more\n\033[0m")

    # FastP don't allow using more than 16 threads
    if int(args.threads) > 16:
        sys.stderr.write("\033[93mWARNING: fastp do not"
                         " allows more than 16 threads."
                         " Threads parameter has been set to 16\n\033[0m")
        args.threads = "16"

    # Check how many threads should been spawned
    CMD.append("--thread")
    CMD.append(args.threads)

    if args.phred64:
        # Specify if Phred +64 scale is being used
        CMD.append("--phred64")

    # Pass the html and json output name
    CMD.append("--html={}.fastp.html".format(args.prefix))
    CMD.append("--json={}.fastp.json".format(args.prefix))

    # Now append parameters defined in config file
    parameters = ReadConfig(args.cfgfile)
    CMD.extend(parameters)

    return CMD


def FastClean(args):
    '''Call fastp'''
    from subprocess import call
    import sys
    import os
    from .History import UpdateHistory

    # Build a list with the command line to call fastP
    fastp_CMD = BuildCMD(args)
    # If verbose is specified, print to the terminal the CMD used to call fastp
    if args.verbose:
        sys.stdout.write("\n{}\n".format(" ".join(fastp_CMD)))

    # When calling fastp, print the output with the summary of the cleaning
    # to PREFIX.cleanlog
    with open("{}.cleanlog".format(args.prefix), "w") as outfh:
        stat = call(fastp_CMD, stdout=outfh, stderr=outfh)
    if stat != 0:
        sys.exit("\033[91mERROR: Pipeline stopped when cleaning with fastp\n"
                 "Stat value: {}\n\033[0m".format(stat))

    # Remove the fastp.html file as this info can be re-generated
    # with the multiqc module
    # and we have also the .json and .cleanlog reports
    os.remove("{}.fastp.html".format(args.prefix))

    UpdateHistory(fastp_CMD, "fastp", args.prefix)

    return 0
