#! /usr/bin/env python3.7
# Copyright (C) 2024 Alvaro Chiner Oms & Miguel Moreno Molina
# adapted from Galo Adrian Goig's script

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

# This is script is thought to get as input a file of type
# "genetic_distances" with following structure
#
#
# Sample1   Distance  Sample2
# ERR1679588    803 G1986
# ERR1679585    2121    ERR1679588
# ERR1679587    2146    ERR1679588
# ERR1679586    2209    ERR1679588
#
# IMPORTANTLY, quotes MUST BE REMOVED in the input file

def ParseDistances(distfile, sep):

    with open(distfile) as infile:
        header = infile.readline()
        if sep == "Space":
            for line in infile:
                line = line.rstrip().split()
                yield line
        elif sep == "Tab":
            for line in infile:
                line = line.rstrip().split("\t")
                yield line


def getClusters(args, threshold):

    clusters = []
    for A, dist, B in ParseDistances(args.distfile, args.sep):
        clustered = False
        dist = int(dist)
        if dist <= threshold:
            if clusters:
                samples_present_in = []
                for i in range(len(clusters)):
                    cluster = clusters[i]
                    if A in cluster:
                        cluster.append(B)
                        samples_present_in.append(i)
                        clustered = True

                    elif B in cluster:
                        clustered = True
                        cluster.append(A)
                        samples_present_in.append(i)

                if not clustered:
                    newcluster = [A, B]
                    clusters.append(newcluster)

                elif len(samples_present_in) > 1:
                    newcluster = []
                    offset = 0
                    for i in samples_present_in:
                        newcluster.extend(clusters[i-offset])
                        del clusters[i-offset]
                        offset += 1
                    clusters.append(newcluster)

            else:
                clusters.append([A, B])
        else:
            if clusters:
                A_present_in = []
                for i in range(len(clusters)):
                    cluster = clusters[i]
                    if A in cluster:
                        A_present_in.append(i)
                        clustered = True

                if not clustered:
                    clusters.append([A])

                elif len(A_present_in) > 1:
                    newcluster = []
                    offset = 0
                    for i in A_present_in:
                        newcluster.extend(clusters[i-offset])
                        del clusters[i-offset]
                        offset += 1
                    clusters.append(newcluster)

                clustered = False
                B_present_in = []
                for i in range(len(clusters)):
                    cluster = clusters[i]
                    if B in cluster:
                        B_present_in.append(i)
                        clustered = True

                if not clustered:
                    clusters.append([B])

                elif len(B_present_in) > 1:
                    newcluster = []
                    offset = 0
                    for i in B_present_in:
                        newcluster.extend(clusters[i-offset])
                        del clusters[i-offset]
                        offset += 1
                    clusters.append(newcluster)

            else:
                clusters.append([A])
                clusters.append([B])

    return clusters


def deDupClusters(clusters):
    dedup = []
    for cluster in clusters:
        cluster = set(cluster)
        cluster = list(cluster)
        dedup.append(cluster)

    return dedup


def writeClusters(clusters, outfile, threshold):
    cont = 1
    file_out = open("{}.clusters_{}.prov".format(outfile, str(threshold)), "w")
    file_out.write("Cluster\tSamples\n")
    for cluster in clusters:
        if len(cluster) > 1: #let's save only clusters with two or more samples
            line = "Cluster_{}\t".format(cont)
            line += "{}".format(cluster[0])
            for sample in cluster[1:]:
                line += ",{}".format(sample)
            line += "\n"
            file_out.write(line)
        cont += 1
    file_out.close()


def writeClustersIrving(clusters, outfile, threshold):
    cont = 1
    file_out = open("{}.clusters_{}.prov".format(outfile, str(threshold)), "w")
    file_out.write("Cluster\tSamples\n")
    for cluster in clusters: 
        if len(cluster) > 1: #let's save only clusters with two or more samples
            line = "Cluster_{}\t".format(cont)
            for sample in cluster:
                file_out.write("{}{}\n".format(line, sample))
        cont += 1
    file_out.close()



def normalizeClusters(infile1, infile2, output):
    import collections

    #we load the first clusters list
    clustersA = {}

    with open(infile1) as infile:
        next(infile)
        for line in infile:
            line = line.strip().split("\t")
            if line[0] not in clustersA.keys():
                clustersA[line[0]] = []
                clustersA[line[0]].append(line[1])
            else:
                clustersA[line[0]].append(line[1])

    #we load the second clusters list
    clustersB = {}

    with open(infile2) as infile:
        next(infile)
        for line in infile:
            line = line.strip().split("\t")
            if line[0] not in clustersB.keys():
                clustersB[line[0]] = []
                clustersB[line[0]].append(line[1])
            else:
                clustersB[line[0]].append(line[1])


    #we now search the old samples in the updated groups
    new_clustersB = {}
    assigned = []
    for cluster, query in clustersA.items():
        for name, samples in clustersB.items():
            if len(set(query)&set(samples)) == len(query):
                new_clustersB[cluster] = samples
                assigned.append(samples)

    #we need to track which are the new clusters at the current threshold
    pending = []
    for group in clustersB.values():
        if group not in assigned:
            pending.append(group)

    current = max([int(x.split("_")[1]) for x in new_clustersB.keys() if "+" not in x])
    for group in pending:
        name = "Cluster_" + str(current + 1)
        current += 1
        new_clustersB[name] = group

    #and finally we merge clusters that have been combined by increasing the threshold
    for cl, p in new_clustersB.copy().items():
        for cl2, p2 in new_clustersB.copy().items():
            if cl != cl2 and p == p2:
                try:
                    new_clustersB[cl+"+"+cl2.split("_")[1]] = new_clustersB.pop(cl)
                    del new_clustersB[cl2]
                except KeyError:
                    pass

    #write result to file
    with open(output, "w") as out:
        out.write("Cluster\tSamples\n")
        for cl, p in new_clustersB.items():
            for m in p:
                out.write(cl + "\t" + m + "\n")


def GetClusters(args):
    '''We will produce 5 files at SNP thresholds 0, 5, 10, 12 and 15.
    Also, we will normalize the clusters naming across thresholds.'''
    import subprocess as sp
    
    
    if args.threshold is None:
        for threshold in [0, 5, 10, 12, 15]:
            print("\033[92mCalculating clusters with threshold {}...\033[00m".format(str(threshold)))

            clusters = getClusters(args, threshold)
            clusters = deDupClusters(clusters)
    
            if not args.output:
                writeClusters(clusters, args.outfile, threshold)
            else:
                writeClustersIrving(clusters, args.outfile, threshold)
                
                #we will only normalize the clusters naming if -osi
        # print("\033[92mNow normalizing clusters...\033[00m")
            
        # normalizeClusters("{}.clusters_0.prov".format(args.outfile),
        #               "{}.clusters_5.prov".format(args.outfile),
        #               "{}.clusters_5.tsv".format(args.outfile))
        # normalizeClusters("{}.clusters_5.tsv".format(args.outfile),
        #               "{}.clusters_10.prov".format(args.outfile),
        #               "{}.clusters_10.tsv".format(args.outfile))
        # normalizeClusters("{}.clusters_10.tsv".format(args.outfile),
        #               "{}.clusters_12.prov".format(args.outfile),
        #               "{}.clusters_12.tsv".format(args.outfile))
        # normalizeClusters("{}.clusters_12.tsv".format(args.outfile),
        #               "{}.clusters_15.prov".format(args.outfile),
        #               "{}.clusters_15.tsv".format(args.outfile))
                      
                #final cleanup
        sp.run("mv {}.clusters_0.prov {}.clusters_0.tsv".format(args.outfile, args.outfile),
                   shell=True, capture_output=True)
          
        #sp.run("rm *prov", shell=True, capture_output=True)
        
    else:
    
        print("\033[92mCalculating clusters with threshold {}...\033[00m".format(str(args.threshold)))
        clusters = getClusters(args, args.threshold)
        clusters = deDupClusters(clusters)
        
        if not args.output:
            writeClusters(clusters, args.outfile, args.threshold)
        else:
            writeClustersIrving(clusters, args.outfile, args.threshold)
            
        sp.run("mv {}.clusters_{}.prov {}.clusters_{}.tsv".format(args.outfile, str(args.threshold), args.outfile, str(args.threshold)),
                   shell=True, capture_output=True)
    
    print("\033[92mDONE!\033[00m")
    
