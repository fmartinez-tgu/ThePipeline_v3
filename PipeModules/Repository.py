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


def Programs():
    ''' Parse .programs_path and return a dictionary with programs paths like
    bwa = /data/Software/bwa/bwa.1.23 '''
    import os

    comm_path = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
    programs = {}
    with open("{}/data/Paths/programs_path".format(comm_path)) as infile:
        for line in infile:
            if not line.startswith("#"):
                line = line.rstrip()
                program, path = line.split("=")
                programs[program] = "{}{}".format(comm_path, path)

    return programs


def Data():
    ''' Parse .data_path and return a dictionary with data files paths like:
    data["reference"] = /data/genome_epi_val/MTB_ancestor.fasta '''
    import os

    comm_path = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
    data = {}
    with open("{}/data/Paths/data_path".format(comm_path)) as infile:
        for line in infile:
            if not line.startswith("#"):
                line = line.rstrip()
                file, path = line.split("=")
                data[file] = path

    return data
