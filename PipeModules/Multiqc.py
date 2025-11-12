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


def Multiqc(args):
    import os
    import subprocess as sp
    from .Repository import Programs
    import sys

    programs = Programs()
    multiqc = programs["multiqc"]

    # multiqc config file
    config = "{}/data/Configs/multiqc.config".format(
            os.path.split(
                os.path.dirname(os.path.abspath(__file__)))[0])

    # to run multiqc wee need to use the specific python interpreter
    # located in MultiQC_env
    cmd_multiqc = [''.join([os.path.dirname(multiqc), '/python']),
                   multiqc, args.folder, '--filename', args.output,
                   '--interactive', '--outdir', '.', '--exclude',
                   'snippy', '-f', '-c', config]

    # run multiqc and show progress
    mem = sp.Popen(cmd_multiqc, stdout=sp.PIPE)
    for line in mem.stdout:
        sys.stdout.write(line.decode())
