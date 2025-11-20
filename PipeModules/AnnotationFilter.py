#! /usr/bin/env python3.7
# Copyrigth (C) 2024 Miguel Moreno Molina

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

def LoadAnnotation(anno_file):
    '''Load the annotation file to store each position to
    be discarded in a list

    This function is meant to load by default
    the H37Rv annotation under
    /data/Databases/MTB_annotation/H37Rv.annotation.tab
    '''
    
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


# def FilterSnps(args):
#     '''Filter SNP files according to the TSV file
#     passed. This TSV file is an annotation file of the
#     H37Rv genome, with the last column indicating if SNPs in 
#     this region must be kept or not'''

#     import os
#     import vcf

#     annotation = LoadAnnotation("{}/data/H37Rv.annotation_new.tsv".format(
#                     os.path.split(
#                         os.path.dirname(
#                             os.path.abspath(__file__)))[0]))

#     with open("{}.EPI.snp.final.annoF".format(args.file_name), "w") as outfile:
#         with open("{}.EPI.snp.final".format(args.prefix)) as infile:
#             header = infile.readline()
#             outfile.write(header)
#             for line in infile:
#                 pos = line.strip().split("\t")[1]
#                 if pos not in annotation:
#                         outfile.write(line)
    
#     vcf_reader = vcf.Reader(open("{}.EPI.snp.vcf".format(args.prefix), 'r'))
#     vcf_output = vcf.Writer(open("{}.EPI.snp.vcf.annoF".format(args.prefix), 'w'), vcf_reader)

#     for record in vcf_reader:
#         if str(record.POS) not in annotation:
#             vcf_output.write_record(record)

#     return 0

def FilterSnps(args):
    '''Filter SNP files of an input file'''

    import os
    import subprocess as sp 
    annotation = LoadAnnotation("/data/ThePipeline_v3/data/H37Rv.annotation_new.tsv".format(
                    os.path.split(
                        os.path.dirname(
                            os.path.abspath(__file__)))[0]))


    with open("{}_annoF".format(args.file_name), "w") as outfile:
        with open("{}".format(args.file_name)) as infile:
            lines_varscan = infile.readlines()
            headers = [x.strip() for x in lines_varscan if "#" in x]
            rest_of_file = [x.strip() for x in lines_varscan if "#" not in x]
            outfile.write("\n".join(headers)+"\n")
            for line in rest_of_file:
                pos = line.strip().split("\t")[1]
                if pos not in annotation:
                        outfile.write(f"{line}\n")