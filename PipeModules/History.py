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

def UpdateHistory(command, program, prefix):
    '''This function will be called to either update a history_tp.txt file where
    with the last command executed, the directory and the date
    '''
    import time
    import os
    from .Version import version
    # Get the date
    date = time.localtime()
    cwd = os.getcwd()
    if program == 'custom':
        v = "ThePipeline3 custom part"
    else:
        v = version(program)

    hist = "{}/{}.history".format(cwd, prefix)
    if os.path.isfile(hist):  # Check wheter history_tp already exists
        if os.stat(hist).st_size != 0:  # Check that is not empty
            with open(hist) as infile:
                history = infile.readlines()  # Read the file
        else:
            history = None
    else:
        history = None

    with open(hist, "w") as outfile:
        if history:
            for line in history:
                outfile.write(line)

        outfile.write("{};".format(program))
        outfile.write(cwd + ";")
        outfile.write("{}-{}-{}@{}:{};".format(
                                                date.tm_mday,
                                                date.tm_mon,
                                                date.tm_year,
                                                date.tm_hour,
                                                date.tm_min
                                                ))
        if program == 'custom':
            outfile.write(command + ";")
        else:
            outfile.write(" ".join(command) + ";")
        outfile.write(v)
        outfile.write("\n")

    return 0
