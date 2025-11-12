#! /usr/bin/env python3.7
# Copyright (C) 2024 Miguel Moreno Molina

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

# Module for determining the MTBC lineage of the samples, based
# on a specific list of lineage definitory SNPs.

BCG_barcode = ["1031219","1101047","1125845","114093","1149361","1171637","1196508","1215076","126042",
               "1281650","128264","1317991","1346228","1354643","1363023","1456926","1489971","1505535",
               "1545125","1606157","1646777","1650880","1672136","1711336","1723697","1749146","1757541",
               "1778259","1798570","1929580","1931093","1995944","2021331","2027617","2095504","2104479",
               "2119208","2120686","2135450","2136619","2144704","2159209","224441","2247905","2254221",
               "2298095","2318692","2326813","2385036","241776","2483513","253032","2543670","2592413",
               "2705340","2729423","2750756","2758803","2782928","2799276","2820456","2888978","2924782",
               "294127","2982422","3033869","3038262","3062237","3076647","3089363","3095342","3138479",
               "3139854","3143158","3180775","3183595","3192053","3207709","324902","3306811","3349019",
               "3423939","3423966","3426021","3431434","3446362","3450822","346034","3489981","3501783",
               "3562107","3726374","3886326","39093","3918328","3958885","3966406","4041134","4048534",
               "4063541","407507","4085376","4090810","4112441","4116617","4117009","4140245","4151717",
               "4167189","423289","4262953","4269689","4330376","4378914","4383446","4409519","484874",
               "494270","505302","50766","534394","558720","559153","562943","579674","580860","587976",
               "643319","65766","659808","664045","682964","69072","763575","781568","822822","84972",
               "855557","8624","864517","87079","872431","916429","946958","961941","983462"]
               
def Typing():
    '''Check if a variant file contains phylogenetic SNPs. If so
    return the lineage or lineages, and the type of infection thereafter'''

    import os, glob
    phylo = {}
    
    markers_file = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0] + "/data/snp_phylo_fixed.tsv"

    with open(markers_file) as infile: #we load the phylogenetic markers file
        next(infile)
        for line in infile:
            line = line.strip().split("\t")
            snpid = line[1]+line[3].replace("/","")
            phylo[snpid] = line[0]
            
    with open("lineage_typing.csv", "w") as outfile:
        outfile.write("Sample,Infection,Typing\n")
        for filename in glob.glob("*.DR.snp.final"):
            positions = []
            with open(filename) as snpfile:
                filename = filename.split("/")[-1]
                try:
                    next(snpfile)
                except:
                    outfile.write(filename.replace(".DR.snp.final","") + ",emptyFile!\n")
                    continue

                lins, markers_freq = [], {}

                for line in snpfile:
                    line = line.strip().split("\t")
                    positions.append(line[1])
                    snpid = line[1] + line[2] + line[6]
                    freq = float(line[4])

                    if snpid in phylo.keys():
                        lins.append(phylo[snpid].replace("lineage","L"))
                        markers_freq[phylo[snpid].replace("lineage","L")] = str(freq).split(".")[0]+"%"

                lins = sorted(lins, key = len)
                #now we determine if we have a clonal or mixed infection
                n_strains = len(set([element.split(".")[0] for element in lins]))   

                if n_strains == 1: #clonal infection
                    lengths = [len(s) for s in lins]
                    deepest = lengths.index(max(lengths))
                    if lins[deepest] == "LBovis":
                        if(all(x in positions for x in BCG_barcode)):
                            lins[deepest] = "LBovis_BCG"

                    exceptions = ["L3.1.1", "L3.1.2", "L1.2.1", "L1.2.2"]
                    
                    if len(lins) > 1: #we check the lineage hierarchy
                        for i in range(len(lins)-1, 0, -1):
                            if lins[i] in exceptions:
                                break
                            if ".".join(lins[i].split(".")[:-1]) != lins[i-1]:
                                lins[deepest] = "WARNING (" + lins[i] + "/" + lins[i-1] + ")"
                    
                    outfile.write(filename.replace(".DR.snp.final","") + ",clonal," + lins[deepest] + "\n")

                elif n_strains == 0: #empty
                    outfile.write(filename.replace(".DR.snp.final","") + ",undetermined,undetermined\n")

                elif n_strains > 1: #mixed infection of different lineages
                    mix = []
                    for strain in set([element.split(".")[0] for element in lins]):
                        this_lin = [x for x in lins if strain in x]
                        lengths = [len(s) for s in this_lin]
                        shallow = lengths.index(min(lengths))
                        mix.append(this_lin[shallow]+"("+markers_freq[this_lin[shallow]]+")")

                    outfile.write(filename.replace(".DR.snp.final","") + ",mixed," + "+".join(mix)+"\n")
                    
