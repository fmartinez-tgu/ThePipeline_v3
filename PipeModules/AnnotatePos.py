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
    
def annotateSingle(snp_id, genes, code, ancestor):

    def revComplement(seq): 
        baseComplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
        letters = list(seq) 
        letters = [baseComplement[base] for base in letters] 
        temp = ''.join(letters)
        return(temp[::-1])


    if snp_id[-1] in ["A","T","C","G"] and snp_id[-2] in ["A","T","C","G"]:
        snp_pos = snp_id[:-2]
        this_gene = ""
        strand = ""
        this_seq = ""
        true_pos = ""

        for gene, coords in genes.items():
            if coords[0] <= int(snp_pos) <= coords[1]:
                this_gene = gene
                strand = coords[2]
                this_seq = ancestor[int(coords[0])-1:int(coords[1])]
                true_pos = int(snp_pos) - coords[0]
        
        if this_gene.startswith("Rv"):
            
            codons = [this_seq[i:i+3] for i in range(0, len(this_seq), 3)]
            this_codon = codons[int(float(true_pos) / 3)]
            codon_pos = int(float(true_pos) % 3)
            
            if strand == "-":
                ncod = str(len(codons) - int((float(true_pos) / 3))).split(".")[0]
                this_codon = revComplement(this_codon)

                mutated_codon = list(this_codon)
                mutated_codon[abs(2-codon_pos)] = revComplement(snp_id[-1])
                mutated_codon = "".join(mutated_codon)

                if code[this_codon] == code[mutated_codon]:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",syn," + this_codon + "," + mutated_codon
                else:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",non-syn," + this_codon + "," + mutated_codon

            else:
                ncod = str((float(true_pos) / 3)+1).split(".")[0]

                mutated_codon = list(this_codon)
                mutated_codon[codon_pos] = snp_id[-1]
                mutated_codon = "".join(mutated_codon)

                if code[this_codon] == code[mutated_codon]:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",syn," + this_codon + "," + mutated_codon
                else:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",non-syn," + this_codon + "," + mutated_codon
        else:
            return this_gene + ",-,intergenic,-,-"
    else:
        return "SNP_ID format error"

def annotateDouble(snp_id, genes, code, ancestor):

    def revComplement(seq): 
        baseComplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
        letters = list(seq) 
        letters = [baseComplement[base] for base in letters] 
        temp = ''.join(letters)
        return(temp[::-1])
        
    if "+" in snp_id:
        snp_pos1 = snp_id.split("+")[0][:-2]
        snp_pos2 = snp_id.split("+")[1][:-2]
        this_gene = ""
        strand = ""
        this_seq = ""
        true_pos1 = ""
        true_pos2 = ""

        for gene, coords in genes.items():
            if coords[0] <= int(snp_pos1) <= coords[1]:
                this_gene = gene
                strand = coords[2]
                this_seq = ancestor[int(coords[0])-1:int(coords[1])]
                true_pos1 = int(snp_pos1) - coords[0]
                true_pos2 = int(snp_pos2) - coords[0]
        
        if this_gene.startswith("Rv"):
            
            codons = [this_seq[i:i+3] for i in range(0, len(this_seq), 3)]
            this_codon = codons[int(float(true_pos1) / 3)]
            codon_pos1 = int(float(true_pos1) % 3)
            codon_pos2 = int(float(true_pos2) % 3)
            
            if strand == "-":
                ncod = str(len(codons) - int((float(true_pos1) / 3))).split(".")[0]
                this_codon = revComplement(this_codon)
            
                mutated_codon = list(this_codon)
                mutated_codon[abs(2-codon_pos1)] = revComplement(snp_id.split("+")[0][-1])
                mutated_codon[abs(2-codon_pos2)] = revComplement(snp_id.split("+")[1][-1])
                mutated_codon = "".join(mutated_codon)

                if code[this_codon] == code[mutated_codon]:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",syn," + this_codon + "," + mutated_codon
                else:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",non-syn," + this_codon + "," + mutated_codon

            else:
                ncod = str((float(true_pos1) / 3)+1).split(".")[0]

                mutated_codon = list(this_codon)
                mutated_codon[codon_pos1] = snp_id.split("+")[0][-1]
                mutated_codon[codon_pos2] = snp_id.split("+")[1][-1]
                mutated_codon = "".join(mutated_codon)

                if code[this_codon] == code[mutated_codon]:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",syn," + this_codon + "," + mutated_codon
                else:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",non-syn," + this_codon + "," + mutated_codon
        else:
            return this_gene + ",-,intergenic"
    else:
        return "SNP_ID format error"
        
        
def annotateTriple(snp_id, genes, code, ancestor):

    def revComplement(seq): 
        baseComplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
        letters = list(seq) 
        letters = [baseComplement[base] for base in letters] 
        temp = ''.join(letters)
        return(temp[::-1])
        
    if "+" in snp_id:
        snp_pos1 = snp_id.split("+")[0][:-2]
        snp_pos2 = snp_id.split("+")[1][:-2]
        snp_pos3 = snp_id.split("+")[2][:-2]
        this_gene = ""
        strand = ""
        this_seq = ""
        true_pos1 = ""
        true_pos2 = ""
        true_pos3 = ""

        for gene, coords in genes.items():
            if coords[0] <= int(snp_pos1) <= coords[1]:
                this_gene = gene
                strand = coords[2]
                this_seq = ancestor[int(coords[0])-1:int(coords[1])]
                true_pos1 = int(snp_pos1) - coords[0]
                true_pos2 = int(snp_pos2) - coords[0]
                true_pos3 = int(snp_pos3) - coords[0]
        
        if this_gene.startswith("Rv"):
            
            codons = [this_seq[i:i+3] for i in range(0, len(this_seq), 3)]
            this_codon = codons[int(float(true_pos1) / 3)]
            codon_pos1 = int(float(true_pos1) % 3)
            codon_pos2 = int(float(true_pos2) % 3)
            codon_pos3 = int(float(true_pos3) % 3)
            
            if strand == "-":
                ncod = str(len(codons) - int((float(true_pos1) / 3))).split(".")[0]
                this_codon = revComplement(this_codon)
            
                mutated_codon = list(this_codon)
                mutated_codon[abs(2-codon_pos1)] = revComplement(snp_id.split("+")[0][-1])
                mutated_codon[abs(2-codon_pos2)] = revComplement(snp_id.split("+")[1][-1])
                mutated_codon[abs(2-codon_pos3)] = revComplement(snp_id.split("+")[2][-1])
                mutated_codon = "".join(mutated_codon)

                if code[this_codon] == code[mutated_codon]:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",syn," + this_codon + "," + mutated_codon
                else:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",non-syn," + this_codon + "," + mutated_codon

            else:
                ncod = str((float(true_pos1) / 3)+1).split(".")[0]

                mutated_codon = list(this_codon)
                mutated_codon[codon_pos1] = snp_id.split("+")[0][-1]
                mutated_codon[codon_pos2] = snp_id.split("+")[1][-1]
                mutated_codon[codon_pos3] = snp_id.split("+")[2][-1]
                mutated_codon = "".join(mutated_codon)

                if code[this_codon] == code[mutated_codon]:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",syn," + this_codon + "," + mutated_codon
                else:
                    return this_gene + "," + code[this_codon] + ncod + code[mutated_codon] + ",non-syn," + this_codon + "," + mutated_codon
        else:
            return this_gene + ",-,intergenic"
    else:
        return "SNP_ID format error"

