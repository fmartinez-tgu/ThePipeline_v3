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

def ReportParser(prefix):
	'''Parse a kraken-report and produce two tables (genus and species) with
    contaminants and its frequency'''
	MTB_species = ["Mycobacterium tuberculosis",
                   "Mycobacterium africanum",
                   "Mycobacterium bovis",
                   "Mycobacterium canettii",
                   "Mycobacterium microti",
                   "Mycobacterium caprae",
                   "Mycobacterium pinnipedii"]
	# Tables are stored in a list so it can be
	# The ordered according to the frequency of
	# contaminants
	genus_table = []
	species_table = []

	# How many reads assigned to genus or species levels
	genus_reads = 0
	species_reads = 0	
	# We cannot assume there is going to be mycobacterium
	# so initialize one variable in case myc_reads and myc_freq
	# is never assigned
	myc_reads = 0
	myc_freq = 0.0
	mtb_reads = 0
	mtb_freq = 0
	
	with  open("{}.kraken.report".format(prefix)) as infile:
		unclass_line = infile.readline()
		unclass_line = unclass_line.rstrip()
		freq, rcov, rassign, level, taxid, taxonomy = unclass_line.split("\t")
		row = [ prefix, "Unclassified", int(rcov), float(freq) ]
		genus_table.append(row)
		species_table.append(row)
		unclass_reads = int(rcov)

		root_line = infile.readline()
		root_line = root_line.rstrip()
		freq, rcov, rassign, level, taxid, taxonomy = root_line.split("\t")
		root_reads = int(rcov)

		total_reads = root_reads + unclass_reads

		for line in infile:
			line = line.rstrip()
			freq, rcov, rassign, level, taxid, taxonomy = line.split("\t")
			taxonomy = taxonomy.strip()  # Remove white spaces
			freq = float(freq)
			rcov = int(rcov)
			if level == "G": # Get unclassified reads
				if taxonomy == "Mycobacterium":
					myc_reads = rcov
					myc_freq = freq
				else:
					row = [ prefix, taxonomy, rcov, freq ]
					genus_table.append(row)
			elif level == "S":
				# We're gonna treat MTB and nonMTB species differently
				# as all MTB species are gonna be classified as MTBC
				if taxonomy not in MTB_species:
					row = [ prefix, taxonomy, rcov, freq ]
					species_table.append(row)
			elif taxonomy == "Mycobacterium tuberculosis complex":
				mtb_reads = rcov
				mtb_freq = freq
				row = [ prefix, taxonomy, rcov, freq ]
				species_table.append(row)
				genus_table.append(row)
			
		#  Take away MTB reads from Mycobacterium genus
	
		myc_reads -= mtb_reads
		myc_freq -= mtb_freq
		myc_row = [ prefix, "Non-TB-Mycobacterium", myc_reads, myc_freq]
		genus_table.append(myc_row)

	# Sort tables by decreasing percentaje 
	genus_table = sorted(genus_table, key=lambda x:x[3], reverse=True)
	species_table = sorted(species_table, key=lambda x:x[3], reverse=True)

	return genus_table, species_table


def AssessContam(prefix):
	'''Write genus and species tables'''
	genus_table, species_table = ReportParser(prefix)

	with open("{}.genus.contaminants".format(prefix), "w") as outfile:
		outfile.write("Sample\tClassification\tReads\tPercentage\n")
		for row in genus_table:
			sample, taxonomy, reads, perc = row
			outfile.write("{}\t{}\t{}\t{}\n".format(sample, taxonomy, reads, perc))

	with open("{}.species.contaminants".format(prefix), "w") as outfile:
		outfile.write("Sample\tClassification\tReads\tPercentage\n")
		for row in species_table:
			sample, taxonomy, reads, perc = row
			outfile.write("{}\t{}\t{}\t{}\n".format(sample, taxonomy, reads, perc))
