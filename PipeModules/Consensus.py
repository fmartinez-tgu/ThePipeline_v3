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
    try:
        vcf_files = VCFtoPandas('positions_total')
    except Exception as e:
        sys.exit("\033[91mERROR: Something went wrong when"
                 " reading the positions_total file. Check all annoF files are correct.\033[0m")
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

    try:
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
    except Exception as e:
        sys.exit("\033[91mERROR: Something went wrong when"
                 " reading the snpeff_concat_deduplicated file. Check all SnpEff files are correct.\033[0m")
    
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

    # Attempt to fix ALT redundancies using snpeff_concat_deduplicated if present
    print("\033[92m\nChecking for redundancies in the SNP table\n\033[00m".format(threads))
    try:
        # call the fixer defined below in this module
        snp_path = "{}.SNP_table.txt".format(outfile)
        snpeff_path = "snpeff_concat_deduplicated"
        if os.path.exists(snpeff_path):
            try:
                fix_snp_table_from_snpeff(snp_path, snpeff_path, out_path=snp_path)
            except Exception as e:
                print(f"Warning: could not auto-fix SNP table: {e}")
    except Exception:
        pass
    
    # Now we read the SNP table again to ensure it's updated with the fixed version if the fixer ran successfully
    #sp.run("rm snpeff_concat*", shell=True)
    SNP_table = pd.read_csv("{}.SNP_table.txt".format(outfile), sep="\t")

    return SNP_table
    
def fix_snp_table_from_snpeff(snp_table_path, snpeff_concat_path, out_path=None):
    """
    Fix the ALT column and annotation-derived columns in a SNP table using the
    SnpEff concatenated deduplicated file.

    - snp_table_path: path to the SNP table (tab-delimited with header 'Position' and 'ALT')
    - snpeff_concat_path: path to the snpeff_concat_deduplicated file used to build the table
    - out_path: optional path where corrected table will be written (if None, overwrite input file)

    The function will:
    - Read the SNP table into a pandas DataFrame
    - Read the SnpEff file into a DataFrame (tab-separated, header)
    - For each POS in the SNP table, gather all ALT alleles present in snpeff file and replace the ALT
      column with a deduplicated, comma-separated sorted ALT list.
    - Populate/replace the Variant_type, Nuc_change and AA_change columns using ANN/INFO fields from SnpEff
      aggregated per (POS, ALT) when possible, or per POS if the specific allele is not present in snpeff.
    - Write corrected table to out_path.
    """
    import pandas as pd
    import os

    if out_path is None:
        out_path = snp_table_path

    # Load SNP table
    snp_df = pd.read_csv(snp_table_path, sep='\t', dtype={
        'Position': 'int64',
        'ALT': 'string'
    })

    # Read snpeff_concat_deduplicated
    if not os.path.exists(snpeff_concat_path):
        raise FileNotFoundError(f"SnpEff file not found: {snpeff_concat_path}")

    snpeff_df = pd.read_csv(snpeff_concat_path, sep='\t', comment='#', header=0, dtype=str)

    # Normalize column names
    snpeff_df.columns = [c.strip() for c in snpeff_df.columns]

    # Ensure POS and ALT columns exist
    if 'POS' not in snpeff_df.columns or 'ALT' not in snpeff_df.columns:
        raise ValueError("snpeff_concat file must contain POS and ALT columns")

    # Build a mapping: Position -> set(ALT)
    pos_to_alts = {}
    pos_to_ref = {}
    # also populate pos_to_ref using the full dataframe (some loops below iterate only subset columns)
    for _, row in snpeff_df.iterrows():
        try:
            p = int(row['POS'])
        except Exception:
            continue
        try:
            ref_val = str(row.get('REF', ''))
            if ref_val and ref_val != '.':
                pos_to_ref.setdefault(p, ref_val)
        except Exception:
            pass

    for _, row in snpeff_df[['POS', 'ALT']].iterrows():
        try:
            p = int(row['POS'])
        except Exception:
            continue
        alt_field = str(row['ALT'])
        # (pos_to_ref already populated above)
        # ALT may contain comma-separated alleles in VCF; split and dedupe
        for a in alt_field.split(','):
            a = a.strip()
            if a == '' or a == '.':
                continue
            pos_to_alts.setdefault(p, set()).add(a)

    # Build annotation aggregates per (POS, ALT) and per POS
    # We'll look for ANN= or the INFO field following SnpEff conventions
    # pos_alt_ann: (pos, alt) -> {'Variant_type': set(), 'Nuc_change': set(), 'AA_change': set(), 'Nuc_to_AA': {nuc: set(aa)}, 'Nuc_to_Vtype': {nuc: set(vtypes)}, 'entries': [(nuc,aa,vtype,gene), ...]}
    pos_alt_ann = {}
    # pos_ann: pos -> {'Variant_type': set(), 'Nuc_change': set(), 'AA_change': set(), 'Nuc_to_AA': {nuc: set(aa)}, 'Nuc_to_Vtype': {nuc: set(vtypes)}, 'entries': [...]}
    pos_ann = {}

    import re

    for _, row in snpeff_df.iterrows():
        try:
            p = int(row['POS'])
        except Exception:
            continue
        alt_field = str(row.get('ALT', ''))
        info_field = str(row.get('INFO', ''))

        # canonical list of ALT values for this VCF row
        alt_values = [a.strip() for a in alt_field.split(',') if a.strip() and a.strip() != '.']

        # Extract ANN=... content if present (SnpEff ANN format)
        ann_content = ''
        if 'ANN=' in info_field:
            try:
                ann_content = info_field.split('ANN=')[1].split(';')[0]
            except Exception:
                ann_content = ''

        # ANN entries are comma-separated; each has pipe-delimited subfields per SnpEff doc
        if ann_content:
            effects = ann_content.split(',')
            # Process each ANN effect and try to assign it to the ALT allele it references
            for eff in effects:
                parts = eff.split('|')
                # parts[0]=allele (usually), parts[1]=impact/effect (Variant type), parts[9]=nuc change, parts[10]=aa change
                raw_allele = parts[0] if len(parts) > 0 else ''
                vval = parts[1] if len(parts) > 1 and parts[1] else ''
                nval = parts[9].replace("c.", "") if len(parts) > 9 and parts[9] else ''
                aval = parts[10].replace("p.", "") if len(parts) > 10 and parts[10] else ''

                # determine allele this ANN effect refers to
                assigned_allele = None
                parsed_mismatch = False
                # Prefer to parse the nucleotide change to extract the ALT base (most reliable)
                if nval:
                    m = re.search(r"[A-Za-z0-9_\-]*?([ACGT])>([ACGT])", nval)
                    if m:
                        alt_base = m.group(2)
                        if alt_base in alt_values:
                            assigned_allele = alt_base
                        else:
                            # we parsed an ALT base but it doesn't match the VCF ALT values -> don't trust ANN allele token
                            parsed_mismatch = True

                # If we couldn't parse from the c. field, and we didn't parse a mismatching ALT, use the allele token from ANN when it matches an ALT
                if assigned_allele is None and not parsed_mismatch and raw_allele and raw_allele != '.':
                    if raw_allele in alt_values:
                        assigned_allele = raw_allele

                # last resort: if the VCF row has only one ALT and we didn't parse a mismatching c. field,
                # assume the effect refers to it. If the c. field parsed but indicated a different ALT,
                # avoid assigning this ANN to the allele.
                if assigned_allele is None and not parsed_mismatch and len(alt_values) == 1:
                    assigned_allele = alt_values[0]

                # Assign effect to allele-specific mapping if we determined an allele
                gene = parts[3] if len(parts) > 3 and parts[3] else ''
                if assigned_allele:
                    key = (p, assigned_allele)
                    entry = pos_alt_ann.setdefault(key, {'Variant_type': set(), 'Nuc_change': set(), 'AA_change': set(), 'Nuc_to_AA': {}, 'Nuc_to_Vtype': {}, 'entries': []})
                    # append ordered entry preserving ANN order (do not duplicate identical tuples)
                    tup = (nval, aval, vval, gene)
                    if tup not in entry['entries']:
                        entry['entries'].append(tup)
                    if vval:
                        entry['Variant_type'].add(vval)
                    if nval:
                        entry['Nuc_change'].add(nval)
                        entry['Nuc_to_AA'].setdefault(nval, set())
                        entry['Nuc_to_Vtype'].setdefault(nval, set())
                        if aval:
                            entry['Nuc_to_AA'][nval].add(aval)
                        if vval:
                            entry['Nuc_to_Vtype'][nval].add(vval)
                    if aval:
                        entry['AA_change'].add(aval)

                # Also aggregate at POS level as a fallback (keep mapping nuc->aa here too)
                entryp = pos_ann.setdefault(p, {'Variant_type': set(), 'Nuc_change': set(), 'AA_change': set(), 'Nuc_to_AA': {}, 'Nuc_to_Vtype': {}, 'entries': []})
                tup_p = (nval, aval, vval, gene)
                if tup_p not in entryp['entries']:
                    entryp['entries'].append(tup_p)
                if vval:
                    entryp['Variant_type'].add(vval)
                if nval:
                    entryp['Nuc_change'].add(nval)
                    entryp['Nuc_to_AA'].setdefault(nval, set())
                    entryp['Nuc_to_Vtype'].setdefault(nval, set())
                    if aval:
                        entryp['Nuc_to_AA'][nval].add(aval)
                    if vval:
                        entryp['Nuc_to_Vtype'][nval].add(vval)
                if aval:
                    entryp['AA_change'].add(aval)

    # Post-process pos_ann to detect positions where all ANN entries report the same nucleotide change.
    # In those cases we will prefer the first ANN entry for fallbacks (first ANN wins for identical changes).
    for p, pdata in list(pos_ann.items()):
        try:
            entries = pdata.get('entries', [])
            # collect all non-empty nuc changes from entries
            nucs = [e[0] for e in entries if e and e[0]]
            if nucs and len(set(nucs)) == 1:
                pdata['entries_all_same_nuc'] = True
                pdata['entries_first'] = entries[0]
            else:
                pdata['entries_all_same_nuc'] = False
        except Exception:
            pdata['entries_all_same_nuc'] = False

    # Now iterate SNP table and fix ALT and annotation columns
    corrected_rows = []
    corrections = []
    for idx, row in snp_df.iterrows():
        pos = int(row['Position'])
        # Prefer the ALT set from snpeff, otherwise keep existing but clean duplicates
        if pos in pos_to_alts:
            alts = sorted(pos_to_alts[pos])
            new_alt = ','.join(alts)
        else:
            # clean the existing ALT field (handle values like 'A,A,C')
            existing = str(row.get('ALT', ''))
            parts = [p.strip() for p in existing.split(',') if p.strip() and p.strip() != '.']
            new_alt = ','.join(sorted(set(parts), key=parts.index))

        # Fill annotation fields: try per (pos, each_alt) and aggregate
        variant_types, nuc_changes, aa_changes = set(), set(), set()
        for a in new_alt.split(','):
            a = a.strip()
            if not a:
                continue
            key = (pos, a)
            if key in pos_alt_ann:
                entry = pos_alt_ann[key]
                variant_types.update(entry['Variant_type'])
                nuc_changes.update(entry['Nuc_change'])
                aa_changes.update(entry['AA_change'])

        # If no per-allele entries found, use POS-level aggregation
        if not variant_types and pos in pos_ann:
            variant_types.update(pos_ann[pos]['Variant_type'])
            nuc_changes.update(pos_ann[pos]['Nuc_change'])
            aa_changes.update(pos_ann[pos]['AA_change'])

        # Get strand early so per-allele matching can consider it
        strand = ''
        if 'Strand' in snp_df.columns:
            try:
                strand = str(snp_df.at[idx, 'Strand'])
            except Exception:
                strand = ''

        # Build per-allele annotation strings aligned with the ALT order
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        per_allele_vtypes = []
        per_allele_nucs = []
        per_allele_aas = []

        # Prepare POS-level ordered lists for fallback allocation (nuc changes)
        pos_nuc_list = []
        pos_v_str = ''
        pos_aa_str = ''
        if pos in pos_ann:
            pos_nuc_list = sorted([x for x in pos_ann[pos]['Nuc_change'] if x])
            pos_v_str = ';'.join(sorted([x for x in pos_ann[pos]['Variant_type'] if x]))
            pos_aa_str = ';'.join(sorted([x for x in pos_ann[pos]['AA_change'] if x]))

        alleles = [a.strip() for a in new_alt.split(',')]
        for i, a in enumerate(alleles):
            if not a:
                per_allele_vtypes.append('')
                per_allele_nucs.append('')
                per_allele_aas.append('')
                continue
            key = (pos, a)
            if key in pos_alt_ann:
                entry = pos_alt_ann[key]
                # Prefer the first ANN entry recorded (ANN order) when available
                if entry.get('entries'):
                    first_entry = entry['entries'][0]
                    n_str = first_entry[0] or ''
                    aa_str = first_entry[1] or ''
                    # prefer vtype from the entry tuple if present
                    if first_entry[2]:
                        v_str = first_entry[2]
                    else:
                        # fall back to nuc->vtype or Variant_type set
                        if n_str and 'Nuc_to_Vtype' in entry and n_str in entry['Nuc_to_Vtype']:
                            v_str = ';'.join(sorted([x for x in entry['Nuc_to_Vtype'][n_str] if x]))
                        else:
                            v_str = ';'.join(sorted([x for x in entry['Variant_type'] if x]))
                else:
                    # Nuc_change per allele: pick a single representative example (first sorted)
                    n_list = sorted([x for x in entry['Nuc_change'] if x])
                    n_str = n_list[0] if n_list else ''
                    aa_str = ';'.join(sorted([x for x in entry['AA_change'] if x]))
                    if n_str and 'Nuc_to_Vtype' in entry and n_str in entry['Nuc_to_Vtype']:
                        v_str = ';'.join(sorted([x for x in entry['Nuc_to_Vtype'][n_str] if x]))
                    else:
                        v_str = ';'.join(sorted([x for x in entry['Variant_type'] if x]))
            else:
                # fallback: prefer POS-level nuc->AA mappings that match the ALT allele
                n_str = ''
                aa_str = ''
                v_str = pos_v_str

                if pos in pos_ann and 'Nuc_to_AA' in pos_ann[pos]:
                    # try to find nuc changes whose ALT base matches this allele
                    matched = False
                    for nuc_change in sorted([x for x in pos_ann[pos]['Nuc_change'] if x]):
                        # parse ALT base from the nuc change (pattern like 907C>T)
                        m_alt = re.search(r">([ACGT])", nuc_change)
                        if not m_alt:
                            continue
                        ann_alt = m_alt.group(1)
                        # For negative strand, ANN c. ALT is in coding strand (complement of genomic ALT)
                        expected = comp.get(a, a) if strand == '-' else a
                        if ann_alt == expected:
                            n_str = nuc_change
                            # If all entries for this POS report the same nuc change, prefer the first ANN entry globally
                            aa_str = ''
                            vtype_choice = ''
                            if pos_ann[pos].get('entries_all_same_nuc') and 'entries_first' in pos_ann[pos]:
                                first = pos_ann[pos]['entries_first']
                                if first and first[0] == nuc_change:
                                    aa_str = first[1] or ''
                                    vtype_choice = first[2] or ''
                            else:
                                # prefer the first pos_ann entry matching this nuc_change to preserve ANN ordering
                                if 'entries' in pos_ann[pos]:
                                    for etup in pos_ann[pos]['entries']:
                                        if etup[0] == nuc_change:
                                            aa_str = etup[1] or ''
                                            vtype_choice = etup[2] or ''
                                            break
                            if not aa_str:
                                aas = pos_ann[pos]['Nuc_to_AA'].get(nuc_change, set())
                                if aas:
                                    aa_str = ';'.join(sorted(aas))
                                else:
                                    aa_str = ';'.join(sorted([x for x in pos_ann[pos]['AA_change'] if x]))
                            # set v_str to the chosen vtype from entries if found, else try Nuc_to_Vtype, else pos_v_str
                            if vtype_choice:
                                v_str = vtype_choice
                            elif n_str and 'Nuc_to_Vtype' in pos_ann[pos] and n_str in pos_ann[pos]['Nuc_to_Vtype']:
                                v_str = ';'.join(sorted([x for x in pos_ann[pos]['Nuc_to_Vtype'][n_str] if x]))
                            else:
                                v_str = pos_v_str
                            matched = True
                            break

                    # If no matching nuc->AA found, fall back to allocating examples sequentially
                    if not matched:
                        if i < len(pos_nuc_list):
                            n_str = pos_nuc_list[i]
                        elif len(pos_nuc_list) > 0:
                            n_str = pos_nuc_list[-1]
                        else:
                            n_str = ''
                        # for AA, attempt to use mapping for that nuc if available
                        aas = pos_ann[pos]['Nuc_to_AA'].get(n_str, set()) if n_str else set()
                        if aas:
                            aa_str = ';'.join(sorted(aas))
                        else:
                            aa_str = ';'.join(sorted([x for x in pos_ann[pos]['AA_change'] if x]))
                else:
                    # previous behavior when no POS-level mapping exists
                    if i < len(pos_nuc_list):
                        n_str = pos_nuc_list[i]
                    elif len(pos_nuc_list) > 0:
                        n_str = pos_nuc_list[-1]
                    else:
                        n_str = ''
                    aa_str = pos_aa_str

            per_allele_vtypes.append(v_str)
            per_allele_nucs.append(n_str)
            per_allele_aas.append(aa_str)

        # final strings: one element per allele, comma-separated; this preserves duplicated AA entries across alleles
        # If the gene is on the negative strand, complement nucleotide letters in the Nuc_change to match coding strand
        strand = ''
        if 'Strand' in snp_df.columns:
            try:
                strand = str(snp_df.at[idx, 'Strand'])
            except Exception:
                strand = ''

        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        if strand == '-':
            def reconstruct_from_ref(n_str, genome_pos, allele):
                """Rebuild nuc change using genomic REF (pos_to_ref) complemented and allele complemented.
                n_str: existing nuc change string (contains cDNA coordinate)
                genome_pos: genomic POS (int)
                allele: ALT allele in genomic orientation (single base)
                """
                if not n_str:
                    return ''
                # try to get cDNA position number from the existing string
                mpos = re.search(r"(\d+)", n_str)
                if not mpos:
                    # fallback to complementing existing pattern
                    def _repl(m):
                        posnum = m.group(1)
                        refb = m.group(2)
                        altb = m.group(3)
                        crefb = comp.get(refb, refb)
                        caltb = comp.get(altb, altb)
                        return f"{posnum}{crefb}>{caltb}"
                    try:
                        return re.sub(r"(\d+)([ACGT])>([ACGT])", _repl, n_str)
                    except Exception:
                        return n_str

                posnum = mpos.group(1)
                # genomic REF base for this POS
                gref = pos_to_ref.get(genome_pos, None)
                # if genomic REF available and allele is a single base, complement both
                if gref and len(gref) == 1 and allele and len(allele) == 1:
                    crefb = comp.get(gref.upper(), gref.upper())
                    caltb = comp.get(allele.upper(), allele.upper())
                    return f"{posnum}{crefb}>{caltb}"
                else:
                    # fallback to complementing the bases found in the existing n_str
                    def _repl2(m):
                        posnum2 = m.group(1)
                        refb = m.group(2)
                        altb = m.group(3)
                        crefb = comp.get(refb, refb)
                        caltb = comp.get(altb, altb)
                        return f"{posnum2}{crefb}>{caltb}"
                    try:
                        return re.sub(r"(\d+)([ACGT])>([ACGT])", _repl2, n_str)
                    except Exception:
                        return n_str

            # rebuild per-allele nuc strings using genomic REF complement and allele complement
            per_allele_nucs = [reconstruct_from_ref(n, pos, alleles[i] if i < len(alleles) else '') if n else '' for i, n in enumerate(per_allele_nucs)]
            # AA changes are reported relative to the protein sequence; keep values but ensure they remain aligned per allele
            # No change to AA strings themselves (SnpEff reports AA on transcript/protein strand)

        # Post-process: ensure Variant_type elements align with the selected nuc_change when possible.
        # If an element equals the aggregated POS-level pos_v_str (e.g. 'missense_variant;synonymous_variant'),
        # prefer the nuc-specific mapping from pos_ann or pos_alt_ann using per_allele_nucs.
        try:
            for i, vt in enumerate(per_allele_vtypes):
                if not vt or (pos in pos_ann and vt == pos_v_str):
                    n_example = per_allele_nucs[i] if i < len(per_allele_nucs) else ''
                    chosen_v = ''
                    if n_example:
                        # check pos-level nuc->vtype
                        if pos in pos_ann and 'Nuc_to_Vtype' in pos_ann[pos] and n_example in pos_ann[pos]['Nuc_to_Vtype']:
                            chosen_v = ';'.join(sorted([x for x in pos_ann[pos]['Nuc_to_Vtype'][n_example] if x]))
                        else:
                            # check allele-specific mapping
                            allele_key = (pos, alleles[i]) if i < len(alleles) else None
                            if allele_key and allele_key in pos_alt_ann and 'Nuc_to_Vtype' in pos_alt_ann[allele_key] and n_example in pos_alt_ann[allele_key]['Nuc_to_Vtype']:
                                chosen_v = ';'.join(sorted([x for x in pos_alt_ann[allele_key]['Nuc_to_Vtype'][n_example] if x]))
                    if chosen_v:
                        per_allele_vtypes[i] = chosen_v
        except Exception:
            # if anything unexpected happens, keep original per_allele_vtypes
            pass

        vtype_str = ','.join(per_allele_vtypes)
        nchange_str = ','.join(per_allele_nucs)
        achange_str = ','.join(per_allele_aas)

        # Update DataFrame row
        old_alt = str(snp_df.at[idx, 'ALT']) if 'ALT' in snp_df.columns else ''
        snp_df.at[idx, 'ALT'] = new_alt
        if old_alt != new_alt:
            corrections.append((pos, old_alt, new_alt))
        if 'Variant_type' in snp_df.columns:
            snp_df.at[idx, 'Variant_type'] = vtype_str
        if 'Nuc_change' in snp_df.columns:
            snp_df.at[idx, 'Nuc_change'] = nchange_str
        if 'AA_change' in snp_df.columns:
            snp_df.at[idx, 'AA_change'] = achange_str

    # Write corrected table
    # Ensure Position_in_resistant_list stays in uppercase 'TRUE' (or empty) to match original output
    try:
        if 'Position_in_resistant_list' in snp_df.columns:
            def _fix_res_val(x):
                # keep empty/NA as empty string
                if pd.isna(x):
                    return ''
                # booleans -> 'TRUE' or ''
                if isinstance(x, bool):
                    return 'TRUE' if x else ''
                xs = str(x).strip()
                if xs == '':
                    return ''
                if xs.upper() == 'TRUE':
                    return 'TRUE'
                if xs.lower() in ('true', 't', '1', 'yes'):
                    return 'TRUE'
                # otherwise keep as-is (string) but uppercase for consistency
                return xs

            snp_df['Position_in_resistant_list'] = snp_df['Position_in_resistant_list'].apply(_fix_res_val)
    except Exception:
        # if anything goes wrong, continue and write the dataframe (don't crash the fixer)
        pass

    snp_df.to_csv(out_path, sep='\t', index=False)

    # Print a brief summary of corrections (positions changed and examples)
    if corrections:
        print(f"Corrected {len(corrections)} position(s)")
        for pos, old, new in corrections[:50]:
            print(f"POS {pos}: {old} => {new}")
    else:
        print("fix_snp_table_from_snpeff: no ALT corrections needed")

    return snp_df

def allFASTAS(table, paths, threads, sample_list):
    '''Generate all consensus fasta of all the samples
    in the analysis folders. Parallelized'''
    import glob
    import multiprocessing as mp
    from functools import partial
    import os

    # set the number of cores to use
    pool = mp.Pool(threads)

    with open("problematic_files.txt", "w") as pf:
        pf.write("List of problematic files during consensus fasta generation:\n")

    if sample_list:
        prefixes = [x for x in paths]
        # pool.map is blocking, so it waits for everything to finish here.
        pool.map(partial(generateFASTA, table), prefixes)
    else:
        all_tasks = [] # Initialize a list to hold all async tasks
        for folder in paths:
            # Only work with samples having annoF vcfs
            annoF_files = glob.glob("{}/*.EPI.snp.final.annoF".format(folder))
            prefixes = [s.replace('.EPI.snp.final.annoF', '') for s in annoF_files]

            # Use extend to add tasks from every folder to the list
            for prefix in prefixes:
                t = pool.apply_async(generateFASTA, args=(table, prefix))
                all_tasks.append(t)
        
        # Now wait for all tasks across all folders to complete
        for task in all_tasks:
            task.get()

    pool.close()
    pool.join()

    with open("problematic_files.txt", "r") as pf:
        lines = pf.readlines()
        if len(lines) == 1:
            os.remove("problematic_files.txt")


def generateFASTA(table, prefix):
    import pandas as pd
    '''Generate the consensus fasta of a sample'''
    from .Calling import VCFtoPandas
    
    # positions to evaluate
    pos = table['Position'].to_list()

    # load wt.txt, .snp.vcf, .indel.vcf, .snp.varscan and .snp.mutect.tab files
    try:
        wt_file = pd.read_csv("{}.wt".format(prefix), sep="\t", header = 0)
    except Exception as e:
        with open("problematic_files.txt", "a") as pf:
            pf.write(f"{prefix}.wt\n")
        print(f"Error loading {prefix}.wt: {e}")
        return
    
    try:
        indel_file = VCFtoPandas("{}.indel.vcf".format(prefix))
    except Exception as e:
        with open("problematic_files.txt", "a") as pf:
            pf.write(f"{prefix}.indel.vcf\n")
        print(f"Error loading {prefix}.indel.vcf: {e}")
        return
    
    try:
        snp_file = pd.read_csv("{}.snp.minos".format(prefix), sep="\t", header=0)
    except Exception as e:
        with open("problematic_files.txt", "a") as pf:
            pf.write(f"{prefix}.snp.minos\n")
        print(f"Error loading {prefix}.snp.minos: {e}")
        return
    
    try:
        lowcov_file = pd.read_csv("{}.lowcov".format(prefix), sep="\t", header=0)
    except Exception as e:
        with open("problematic_files.txt", "a") as pf:
            pf.write(f"{prefix}.lowcov\n")
        print(f"Error loading {prefix}.lowcov: {e}")
        return
    
    try:
        varscan_file = pd.read_csv("{}.snp.varscan".format(prefix), sep="\t", header=0)
    except Exception as e:
        with open("problematic_files.txt", "a") as pf:
            pf.write(f"{prefix}.snp.varscan\n")
        print(f"Error loading {prefix}.snp.varscan: {e}")
        return
    
    try:
        mutect_file = pd.read_csv("{}.snp.mutect.tab".format(prefix), sep="\t", header=0)
    except Exception as e:
        with open("problematic_files.txt", "a") as pf:
            pf.write(f"{prefix}.snp.mutect.tab\n")
        print(f"Error loading {prefix}.snp.mutect.tab: {e}")
        return


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

    """Write sequence as a single line to avoid concatenation issues."""
    # Ensure seq is a string and remove any pre-existing newlines/whitespace
    clean_seq = str(seq).replace('\n', '').replace('\r', '').strip()

    filename = f"{id}.fas"
    with open(filename, 'w') as f:
        f.write(f">{os.path.basename(id)}\n")
        f.write(f"{clean_seq}\n")


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

    # First, we check whether all the fastas have the same length. If not, we move the problematic fastas to a separate folder
    # and proceed with the remaining fastas to generate the multifasta
    
    fasta_length = table.shape[0] # number of positions in the SNP table, which equals the desired length of each fasta

    lengths = {}
    if sample_list:
        for file in paths:
            with open("{}.fas".format(file), "r") as f:
                lines = f.readlines()
                seq = ''.join([line.strip() for line in lines if not line.startswith('>')])
                lengths[file] = len(seq)
    else:
        for folder in paths:
            fasta_files = sp.run("ls {}/*.fas".format(folder),
                                shell=True, capture_output=True).stdout.decode().strip().split('\n')
            for file in fasta_files:
                with open(file, "r") as f:
                    lines = f.readlines()
                    seq = ''.join([line.strip() for line in lines if not line.startswith('>')])
                    lengths[file] = len(seq)
    
    # Check if all lengths are the same. If not, move the problematic files to a separate folder
    unique_lengths = set(lengths.values())
    if len(unique_lengths) > 1:
        print("\033[93mWARNING: Not all fastas have the same length. Problematic fastas:\033[00m")
        for file, length in lengths.items():
            if length != fasta_length:
                print(f"\033[93m{file} (length: {length})\033[00m")
                sp.run("mv {} ztemp_individual_fastas".format(file), shell=True, capture_output=True)
        print("\033[93mThe problematic fastas have been moved to the ztemp_individual_fastas folder.\033[00m")
        print("\033[93mProceeding to generate multifasta with the remaining fastas.\033[00m")


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

    # Calculate the invariants but considering a reference size without the masked regions
    A = 724576 - table_snpsites_nr[table_snpsites_nr['WT'] == 'A'].shape[0]
    G = 1367600 - table_snpsites_nr[table_snpsites_nr['WT'] == 'G'].shape[0]
    T = 725052 - table_snpsites_nr[table_snpsites_nr['WT'] == 'T'].shape[0]
    C = 1373340 - table_snpsites_nr[table_snpsites_nr['WT'] == 'C'].shape[0]

    # write to file
    with open('{}_invariants_annoF.txt'.format(outfile), 'w') as f:
        f.write("##Number of non-masked invariant sites for "
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

    # Folder to store individual fasta files given any error during the script. If the multifasta is created successfully,
    # this folder is removed. Otherwise, it remains for debugging purposes.
    sp.run("mkdir ztemp_individual_fastas", shell=True, capture_output=True)

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
            sp.run("mv {}/*.fas ztemp_individual_fastas".format(folder),
                shell=True, capture_output=True)

        if args.sample_list:
            prefixes = [x for x in paths]
            for prefix_to_remove in prefixes:
                sp.run("mv {}.fas ztemp_individual_fastas".format(prefix_to_remove),shell=True,capture_output=True)

        sys.exit("\033[91mERROR: Something went wrong when"
                 " generating the individual FASTAs"
                 " in parallel.\033[0m")

    # Generates several multifastas and information associated
    try:
        print("\033[92m\nGenerating multifasta files\n\033[00m")
        multifastas(table, args.outfile, snpsites, paths, args.sample_list)

    except:
        for folder in paths:
            sp.run("mv {}/*.fas ztemp_individual_fastas".format(folder),
                shell=True, capture_output=True)

        if args.sample_list:
            prefixes = [x for x in paths]
            for prefix_to_remove in prefixes:
                sp.run("mv {}.fas ztemp_individual_fastas".format(prefix_to_remove),shell=True,capture_output=True)
        sys.exit("\033[91mERROR: Something went wrong when"
                 " performing the snp-sites step.\033[0m")


    # Remove individual fasta files

    for folder in paths:
        sp.run("rm {}/*.fas".format(folder),
            shell=True, capture_output=True)
    
    if args.sample_list:
        prefixes = [x for x in paths]
        for prefix_to_remove in prefixes:
            sp.run("rm {}.fas".format(prefix_to_remove),shell=True,capture_output=True)
    
    if not os.listdir("ztemp_individual_fastas"):
    # The list of contents is empty, so we can remove it
        sp.run("rm -r ztemp_individual_fastas", shell=True, capture_output=True)


