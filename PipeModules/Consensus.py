#! /usr/bin/env python3.7
# Copyright (C) 2025 Alvaro Chiner-Oms & Miguel Moreno Molina & Francisco Jose Martínez 

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

# Module for generating consensus multiFASTA

def annotate_position(pos, snpeff_df, vcf_files, Annotation, catalog):
    import pandas as pd

    result = {
        'Position': pos,
        'WT': '',
        'ALT': '',
        'Gene_start': '',
        'Gene_end': '',
        'Strand': '',
        'Gene_length': '',
        'Synonym': '',
        'Rv_number': '',
        'DR_gene': '',
        'Sanger_annot': '',
        'Essentiality': '',
        'Annotation': '',
        'Variant_type': '',
        'Nuc_change': '',
        'AA_change': '',
        'Position_in_resistant_list': ''
    }

    try:
        row = snpeff_df[snpeff_df['POS'] == pos]
        alt = ','.join(row['ALT'].drop_duplicates().tolist())
        ref = row['REF'].drop_duplicates().tolist()[0]

        vtype = []
        n_change = []
        a_change = []

        for _, r in row.iterrows():
            fields = r['INFO'].split('|')
            if len(fields) > 10:
                vtype.append(fields[1])
                n_change.append(fields[9])
                a_change.append(fields[10])

        vtype = ','.join(vtype)
        n_change = ','.join(n_change)
        a_change = ','.join(a_change)

        Rv = row['INFO'].tolist()[0].split('|')[3]
        overlap = ''
        anno_row = Annotation[(Annotation['Start'] <= pos) & (Annotation['End'] >= pos)]

        if anno_row.shape[0] > 1:
            anno_row = anno_row[anno_row['Synonym'] == Rv]
            overlap = '+OL'

    except:
        anno_row = Annotation[(Annotation['Start'] <= pos) & (Annotation['End'] >= pos)]
        overlap = '' 
        if anno_row.shape[0] > 1:
            Rv = anno_row['Synonym'].iloc[0]
            anno_row = anno_row[anno_row['Synonym'] == Rv]
            overlap = '+OL'
        else:
            Rv = anno_row['Synonym'].iloc[0]

        ref = vcf_files[vcf_files['POS'] == pos]['REF'].drop_duplicates().tolist()[0]
        alt = ','.join(vcf_files[vcf_files['POS'] == pos]['ALT'].drop_duplicates().tolist())
        vtype = ''
        n_change = ''
        a_change = ''

    if 'IG' in Rv:
        vtype = ''
        n_change = ''
        a_change = ''

    DR_gene = ''
    try:
        dr_field = anno_row['Possible.DR.phage.PEPPE.'].iloc[0]
        if isinstance(dr_field, str) and 'DR-' in dr_field:
            DR_gene = dr_field
    except:
        pass

    res_pos = 'TRUE' if pos in catalog['Position'].tolist() else ''

    result.update({
        'Position': pos,
        'WT': ref,
        'ALT': alt,
        'Gene_start': anno_row['Start'].iloc[0],
        'Gene_end': anno_row['End'].iloc[0],
        'Strand': anno_row['Strand'].iloc[0],
        'Gene_length': anno_row['Length'].iloc[0],
        'Synonym': anno_row['Gene'].iloc[0],
        'Rv_number': Rv + overlap,
        'DR_gene': DR_gene,
        'Sanger_annot': anno_row['sanger_anot'].iloc[0],
        'Essentiality': anno_row['essential.'].iloc[0],
        'Annotation': anno_row['Annotation'].iloc[0],
        'Variant_type': vtype,
        'Nuc_change': n_change.replace('c.', ''),
        'AA_change': a_change.replace('p.', ''),
        'Position_in_resistant_list': res_pos
    })

    return result

def generateSNPtable(paths, outfile, sample_list, threads):
    '''Generate a non redundant list of SNPs, with all the samples of interest
    The SNPs include information about the gene and are annotated'''
    import sys
    import glob
    import pandas as pd
    import io
    import os
    from .Calling import VCFtoPandas
    import subprocess as sp

    SNP_table = pd.DataFrame(columns=['Position', 'WT', 'ALT',
                                          'Gene_start', 'Gene_end', 'Strand',
                                          'Gene_length', 'Synonym',
                                          'Rv_number', 'DR_gene',
                                          'Sanger_annot', 'Essentiality',
                                          'Annotation', 'Variant_type',
                                          'Nuc_change', 'AA_change',
                                          'Position_in_resistant_list'])
    # Load annotation file
    with open("/data/ThePipeline_v3/data/H37Rv.annotation_new.tsv".format(
              os.path.split(
                os.path.dirname(
                    os.path.abspath(__file__)))[0]),
              'r') as f:
        lines = [li for li in f]
    Annotation = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

    # Load catalog of DR mutations
    catalog = pd.read_csv(
                "/data/ThePipeline_v3/data/resistance_positions.csv".format(
                    os.path.split(
                        os.path.dirname(
                            os.path.abspath(__file__)))[0]))

    # Read all VCF.annoF files and concatenate all in positions_total
    vcf_files = pd.DataFrame()
    try:
        sp.run("rm positions_total",
               shell=True)
    except:
        pass

    sp.run('printf "#CHROM\tPOS\tREF\tALT\tVarFreq\tCov_total\tVarAllele\n" >> positions_total',
           shell=True)

    if sample_list:
        files = [x+".EPI.snp.final.annoF" for x in paths]
        # show error if no annoF files
        if len(files) < 1:
            sys.exit("\033[91mERROR: No .annoF "
                     "files in {} folder\033[0m".format(folder))

        for file in files:
            sp.run("cat {} | grep -v '#' >>"
                " positions_total".format(file),
                shell=True)
    else:
        for folder in paths:
            files = glob.glob("{}/*EPI.snp.final.annoF".format(folder))

            # show error if no annoF files
            if len(files) < 1:
                sys.exit("\033[91mERROR: No .annoF "
                        "files in {} folder\033[0m".format(folder))
            else:
                for file in files:
                    sp.run("cat {} | grep -v '#' >> positions_total".format(file), shell=True)

    
    # Later, join in a single pandas dataframe
    vcf_files = VCFtoPandas('positions_total')
    vcf_files.sort_values('POS', inplace=True)
    vcf_files.drop_duplicates(subset=["POS", "ALT"], inplace=True)
    sp.run("rm positions_total", shell=True)


    print("\033[92m\nConcatenating SnpEff files\n\033[00m")
    sp.run("rm snpeff_concat*", shell=True)
    sp.run('printf "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n" >> snpeff_concat',
           shell=True)

    if sample_list:
        snpeff_files = [x+".final_sin_wt_complemented_annotSnpEff.vcf" for x in paths]
        for file in snpeff_files:
            sp.run("cat {} | grep -v '#' >>"
                   " snpeff_concat".format(file), shell=True)

    else:
        for folder in paths:
            snpeff_files = glob.glob("{}/*.final_sin_wt_complemented_annotSnpEff.vcf".format(folder))
            for file in snpeff_files: 
                sp.run("cat {} | grep -v '#' >>"
                   " snpeff_concat".format(file), shell=True)
    
    print("\033[92m\nSorting and deduplicating SnpEff files using {} threads\n\033[00m".format(threads))
    sp.run("sort --parallel={} -k2,2n -k5,5 snpeff_concat | awk '!seen[$2,$5]++' >>"
           " snpeff_concat_deduplicated".format(threads), shell=True)

    snpeff_df = pd.read_csv("snpeff_concat_deduplicated", sep="\t", comment='#', header=0, dtype={
    "CHROM": "string",
    "POS": "int32",
    "ID": "string",
    "REF": "string",
    "ALT": "string",
    "QUAL": "string",
    "FILTER": "string",
    "INFO": "string",
    "FORMAT": "string",
    "SAMPLE": "string"
    })
    
    
    import multiprocessing as mp
    from functools import partial

    print("\033[92m\nAnnotating SNP_table using {} threads\n\033[00m".format(threads))

    # Prepare arguments for parallelization
    annotate_partial = partial(annotate_position,
                               snpeff_df=snpeff_df,
                               vcf_files=vcf_files,
                               Annotation=Annotation,
                               catalog=catalog)

    with mp.Pool(processes=threads) as pool:
        results = pool.map(annotate_partial, vcf_files['POS'].unique())

    SNP_table = pd.DataFrame(results)
    SNP_table.to_csv("{}.SNP_table.txt".format(outfile), sep="\t", index=False)
    SNP_table.index = range(0, SNP_table.shape[0], 1)

    sp.run("rm snpeff_concat*", shell=True)

    return SNP_table



def allFASTAS(table, paths, threads, sample_list):
    '''Generate all consensus fasta of all the samples
    in the analysis folders. Paralelized'''
    import glob
    import multiprocessing as mp

    # set the number of cores to use
    pool = mp.Pool(threads)

    if sample_list:
        prefixes = [x for x in paths]
        pool.map(partial(generateFASTA, table), prefixes)

    else:
        for folder in paths:
            # Only work with samples having annoF vcfs, as were the ones
            # used to generate de SNP table
            annoF_files = glob.glob("{}/*.EPI.snp.final.annoF".format(folder))
            prefixes = [s.replace('.EPI.snp.final.annoF', '') for s in annoF_files]

            # run in parallel for each folder
            tasks = [pool.apply_async(generateFASTA,
                            args=(table, prefix)) for prefix in prefixes]
    for task in tasks:
        task.get()
    pool.close()
    pool.join()


def generateFASTA(table, prefix):
    import pandas as pd
    '''Generate the consensus fasta of a sample'''
    from .Calling import VCFtoPandas
    
    # positions to evaluate
    pos = table['Position'].to_list()

    # load wt.txt, .snp.vcf, .indel.vcf, .snp.varscan and .snp.mutect.tab files
    wt_file = pd.read_csv("{}.wt".format(prefix), sep="\t")
    indel_file = VCFtoPandas("{}.indel.vcf".format(prefix))
    snp_file = pd.read_csv("{}.snp.minos".format(prefix), sep="\t", header=0)
    lowcov_file = pd.read_csv("{}.lowcov".format(prefix), sep="\t", header=0)
    varscan_file = pd.read_csv("{}.snp.varscan".format(prefix), sep = "\t", header=0)
    mutect_file = VCFtoPandas("{}.snp.mutect.tab".format(prefix))

    # keep only deletions (>1 in REF)
    # and add a new colum with the range of the deletion
    # in .indel.vcf
    indel_file = indel_file[indel_file['REF'].apply(len) > 1]
    indel_file['TO'] = indel_file['POS'] + indel_file['REF'].apply(len)
    
    # Convert indel ranges to a list of tuples for faster checking
    indel_ranges = list(zip(indel_file['POS'], indel_file['TO']))
    
    # Create a dictionary for each file for faster lookups
    wt_set = set(wt_file['Position'].to_list())
    lowcov_set = set(lowcov_file['Position'].to_list())

    # For Minos, VarScan and Mutect2 files, we'll remove those variants with a frequency lower than 90% to avoid rescuing non-fixed positions 
    varscan_file = varscan_file[pd.to_numeric(varscan_file['VarFreq'].str.replace('%', ''), errors='coerce') >= 90]
    mutect_file = mutect_file[pd.to_numeric(mutect_file['VarFreq'], errors='coerce') >= 90]

    # After filtering, we create the dictionary 
    snp_dict = snp_file.set_index('Position')['Cons'].to_dict()
    varscan_dict = varscan_file.set_index('Position')['Cons'].to_dict()
    mutect_dict = mutect_file.set_index('Position')['Cons'].to_dict()
    
    # We generate a WT sequence for the sample whose positions will be modified if necessary in the following for loop
    fasta_seq = table['WT'].to_list()

    def is_in_indel_ranges(position, ranges):
        for start, end in ranges:
            if start < position < end:
                return True
        return False
    
    # For each position in the SNP table, check if it appears in in any of the files.
    # In that case, keep the info from that file. 

    for i, p in enumerate(pos):
        if p in snp_dict:
            fasta_seq[i] = snp_dict[p]
        elif p in varscan_dict:
            fasta_seq[i] = varscan_dict[p]
        elif p in mutect_dict:
            fasta_seq[i] = mutect_dict[p]    
        elif p in wt_set:
            continue
        elif is_in_indel_ranges(p, indel_ranges):
            fasta_seq[i] = "-"
        elif p in lowcov_set:
            fasta_seq[i] = "N"
        else:
            fasta_seq[i] = "?"

    # write FASTA
    write_fasta(''.join(fasta_seq), prefix)


def write_fasta(seq, id, wrap=80):
    """Write sequences to a fasta file.

    Parameters
    ----------
    seq : str
        nucleotide sequence
    id : str
        Prefix for naming the fasta file.
    wrap: int
        Number of AA/NT before the line is wrapped.
    """
    import os

    with open('{}.fas'.format(id), 'w') as f:
        f.write('>{}\n'.format(os.path.basename(id)))
        for i in range(0, len(seq), wrap):
            f.write('{}\n'.format(seq[i:i + wrap]))
    f.close()


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    https://www.biostars.org/p/710/
    """
    from itertools import groupby

    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)


def multifastas(table, outfile, snpsites, paths, sample_list):
    '''Module for performing some useful operations
    over the multifasta'''
    import subprocess as sp
    import pandas
    import io

    # First, generate the multifasta file

    if sample_list:
        for file in paths:
           sp.run("cat {}*.fas >> {}.mf.fasta".format(file,outfile),
            shell=True, capture_output=True)
    else:
        for folder in paths:
            sp.run("cat {}/*.fas >> {}.mf.fasta".format(folder,outfile),
                shell=True, capture_output=True)
        
    # Later the gapped alignment, i.e. we'll replace non-ACGT by gaps
    sp.run("cat {}.mf.fasta | sed '/^>/! s/[^ACGTacgt]/-/g' >"
           " {}.mf_gap.fasta".format(outfile, outfile),
           shell=True, capture_output=True)
    
    # run snp-sites and obtain the associated FASTA
    cmd_fastasnpsites = [snpsites,
                         "-o",
                         "{}.mf_gap.snp-sites.fasta".format(outfile),
                         "{}.mf_gap.fasta".format(outfile)]
    sp.call(cmd_fastasnpsites)
    
    # run snp-sites and obtain the associated VCF
    sp.run("snp-sites -v {0}.mf_gap.fasta > {1}_temp.mf_gap.SNP_table.txt".format(outfile, outfile), shell=True, capture_output=True)  
    
    with open("{0}_temp.mf_gap.SNP_table.txt".format(outfile), "r+") as input_temp_vcf:
        lines_temp_vcf = input_temp_vcf.readlines()
        variants = [line for line in lines_temp_vcf if "#" not in line]
        variants_list = []
        for variant in variants:
            variants_list.append(variant.split("\t")[1])
 
    with open("{0}.SNP_table.txt".format(outfile), "r+") as input_original_SNP_table:
        lines_SNP_table = input_original_SNP_table.readlines()
        lines_to_be_kept = []
        header = lines_SNP_table[0]
        for variant in variants_list:
            lines_to_be_kept.append(lines_SNP_table[int(variant)])
   
    with open("{0}.SNP_table.snp-sites.txt".format(outfile), "w+") as output_SNP_table:
        output_SNP_table.write(header)
        for line_add in lines_to_be_kept:
            output_SNP_table.write(line_add)
        
    table_snpsites = pandas.read_csv("{0}.SNP_table.snp-sites.txt".format(outfile), sep='\t', dtype={'Position_in_resistant_list': str})

    res_pos = table_snpsites[table_snpsites['Position_in_resistant_list'] == 'TRUE'].index.tolist()

    # read the snp-sites FASTA
    fasta_snpsites = fasta_iter("{}.mf_gap.snp-sites.fasta".format(outfile))

    # create new fasta
    f = open('{}.mf_gap.snp-sites.no-resis.fasta'.format(outfile), 'w')

    # remove resis positions and write fasta
    for header, seq in fasta_snpsites:
        count = 0
        for remove_pos in res_pos:
            remove_pos = remove_pos - count
            seq = seq[:remove_pos] + seq[remove_pos+1:]
            count += 1
        f.write('>{}\n'.format(header))
        f.write('{}\n'.format(seq))
    f.close()

# finally, generates a new SNP_table, without the removed positions
    table_snpsites_nr = table_snpsites.drop(res_pos, axis=0)
    table_snpsites_nr.index = range(0, table_snpsites_nr.shape[0], 1)
    table_snpsites_nr.to_csv("{}.SNP_table.snp-sites.no-resis.txt".format(
                             outfile),
                             sep="\t", index=False)

    # calculate the invariants
    A = 758359 - table_snpsites_nr[table_snpsites_nr['WT'] == 'A'].shape[0]
    G = 1444793 - table_snpsites_nr[table_snpsites_nr['WT'] == 'G'].shape[0]
    T = 758224 - table_snpsites_nr[table_snpsites_nr['WT'] == 'T'].shape[0]
    C = 1450156 - table_snpsites_nr[table_snpsites_nr['WT'] == 'C'].shape[0]

    # write to file
    with open('{}_invariants.txt'.format(outfile), 'w') as f:
        f.write("##Number of invariant sites for "
                "{}.mf_gap.snp-sites.no-resis.fasta\n".format(outfile))
        f.write("##A,C,G,T\n")
        f.write("{},{},{},{}\n".format(A, C, G, T))
    f.close()
    
    # Remove the temporary VCF file made with snp-sites
    sp.run("rm {0}_temp.mf_gap.SNP_table.txt".format(outfile), shell=True, capture_output=True)

def Consensus(args):
    import os
    import sys
    from .Repository import Programs
    import datetime
    import subprocess as sp

    paths = []
    programs = Programs()

    # if not prefix defined, use current date
    if not args.outfile:
        e = datetime.datetime.now()
        args.outfile = e.strftime("%Y%m%d")

    snpsites = "/data/ThePipeline_v3/Programs/snp-sites/src/snp-sites"

    # Check if the input is a list or a path
    if args.sample_list:
        with open(args.paths[0], "r+") as input_file:
            for line in input_file:
                paths.append(line.strip())

        try:
            print("\033[92m\nCreating SNP table\n\033[00m")
            table = generateSNPtable(paths, args.outfile, args.sample_list, args.threads)
            print("Non redundant SNP list with"
                " associated information saved in"
                " {}.SNP_table.txt".format(args.outfile))
        except Exception as e:
            print(e)
            sys.exit("\033[91mERROR: Something failed "
                    "when generating the SNP table."
                    " Check that .annoF files are present in the"
                    " folder and that all of them are annotated.\033[0m")
    else:    
    # Check that paths are correct
        for f in args.paths:
            if os.path.exists(os.path.abspath(f)):
                paths.append(os.path.abspath(f))
            else:
                sys.exit("\033[91mERROR:{} path does not exists.\033[0m".format(f))
        try:
            print("\033[92m\nCreating SNP table\n\033[00m")
            table = generateSNPtable(paths, args.outfile, args.sample_list, args.threads)
            print("Non redundant SNP list with"
                " associated information saved in"
                " {}.SNP_table.txt".format(args.outfile))
        except Exception as e:
            print(e)
            sys.exit("\033[91mERROR: Something failed "
                    "when generating the SNP table."
                    " Check that .annoF files are present in the"
                    " folder and that all of them are annotated.\033[0m")

    # Generates one consensus for each sample
    try:
        print("\033[92m\nGenerating individual fastas\n\033[00m")
        allFASTAS(table, paths, args.threads, args.sample_list)
        print("\033[92m\nIndividual fasta files generated\n\033[00m")

    except Exception as e:
        print(e)
        for folder in paths:
            sp.run("rm {}/*.fas".format(folder),
                shell=True, capture_output=True)

        if args.sample_list:
            prefixes = [x for x in paths]
            for prefix_to_remove in prefixes:
                sp.run("rm {}.fas".format(prefix_to_remove),shell=True,capture_output=True)

        sys.exit("\033[91mERROR: Something went wrong when"
                 " generating the individual FASTAs"
                 " in parallel.\033[0m")

    # Generates several multifastas and information associated
    try:
        print("\033[92m\nGenerating multifasta files\n\033[00m")
        multifastas(table, args.outfile, snpsites, paths, args.sample_list)

    except:
        for folder in paths:
            sp.run("rm {}/*.fas".format(folder),
                shell=True, capture_output=True)

        if args.sample_list:
            prefixes = [x for x in paths]
            for prefix_to_remove in prefixes:
                sp.run("rm {}.fas".format(prefix_to_remove),shell=True,capture_output=True)
        sys.exit("\033[91mERROR: Something went wrong when"
                 " performing the snp-sites step.\033[0m")

    sp.run("rm *.fas", shell=True, capture_output=True)

    # Remove individual fasta files

    for folder in paths:
        sp.run("rm {}/*.fas".format(folder),
            shell=True, capture_output=True)
    
    if args.sample_list:
        prefixes = [x for x in paths]
        for prefix_to_remove in prefixes:
            sp.run("rm {}.fas".format(prefix_to_remove),shell=True,capture_output=True)

