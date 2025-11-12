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

# Module for calculating the pairwise genetic distances (in number of SNPs)
# between eact pair of sequences in a FASTA file.

import pandas


def parse_output_matrix(indexes):
    return (snames[indexes[0]],
            int(d[indexes[0], indexes[1]]),
            snames[indexes[1]])


def callback_function(res):
    global result
    result = pandas.DataFrame(res, columns=['Sequence_1',
                                            'Distance',
                                            'Sequence_2'])


def Distances(args):
    import sys
    import datetime
    import itertools
    import multiprocessing as mp

    global snames, d, result
    # check that pairsnp is installed
    try:
        from pairsnp import calculate_snp_matrix, calculate_distance_matrix
    except ModuleNotFoundError:
        sys.exit("\033[91mPIPELINE ERROR: Module pairsnp needed!"
                 " Install pairsnp from"
                 " https://github.com/gtonkinhill/pairsnp\033[0m")

    # if not prefix defined, use current date
    if not args.outfile:
        e = datetime.datetime.now()
        args.outfile = e.strftime("%Y%m%d")

    # calculate distances pairsnp python
    sparse_matrix, consensus, snames = calculate_snp_matrix(args.fasta)
    d = calculate_distance_matrix(sparse_matrix, consensus, "dist", False)

    print("Distance matrix calculated\n")
    print("Now formatting output...\n")

    # format the output
    pool = mp.Pool(args.threads)
    params = list(itertools.combinations(range(d.shape[0]), 2))
    pool.map_async(parse_output_matrix, params,
                   callback=callback_function)

    pool.close()
    pool.join()

    result = result.sort_values(by=['Distance'])
    result.to_csv("{}.genetic_distances.tsv".format(args.outfile), sep="\t",
                  index=False)
    print("Results printed to {}.genetic_distances.tsv"
          " file. \n".format(args.outfile))

    if args.limit >= 0:
        result = result[result['Distance'] <= args.limit]
        result.to_csv("{}.genetic_distances_{}.tsv".format(
                                args.outfile, args.limit),
                      sep="\t",
                      index=False)
        print("Results printed to {}.genetic_distances_{}.tsv"
              " file. \n".format(args.outfile, args.limit))
