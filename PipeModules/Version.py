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

# Este modulo es para obtener la version de los programas. He intentado hacerlo
# de un modo fancy, cogiendo la version del output de cada programa, pero cada
# uno es de su padre y su madre asa que habra que mantener la lista actualizada


def version(program):
    '''Read the program versions from data/Configs/sofware_versions.txt'''
    import os

    path = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]

    versions = {}
    with open("{}/data/Configs/software_versions.txt".format(path)) as infile:
        for line in infile:
            soft, version = line.rstrip().split("=")
            versions[soft] = version

    return versions[program]


def ThePipelineV():
    '''This is a suuperspecial module that stores major Pipeline versions. It
    must return both the name of the version and a list of major changes.
    The changes of consecutive versions will
    be referred to the original version
    '''

    TP_version = "ThePipeline3"
