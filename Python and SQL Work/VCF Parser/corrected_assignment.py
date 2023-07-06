"""
Open a GFF/VCF/FASTA file and checks the QUAl attribute of each sample
"""
# IMPORTS
import vcf
import gffutils
import sqlite3
import seaborn as sns
import matplotlib.pyplot as plot
import pandas as pd
from Bio.Seq import Seq
import Bio
import math
import logging
import sys
import os
import warnings
import argparse

# Log Initialiser
logging.getLogger()
logging.basicConfig(filename='2117832_Log.txt', level=logging.INFO, filemode='w',
                    format='%(asctime)s | %(message)s', datefmt='%H:%M:%S')

# Command line interface
interface = argparse.ArgumentParser(prog='VCF file parser',
                                    description='This program takes a VCF, GFF, and fasta file to generate'
                                                'a polished final table of the SNPs stored in the VCF file'
                                                'and the effect of those mutations. A bar plot for the counts'
                                                'of non-coding, non-synonymous, and synonymous mutations.')

interface.add_argument('filename', nargs='*', help='Input files can only be either GFF/VCF/fasta formats')

interface_args = interface.parse_args()

# Input File handler
file_list_str = ''
file_list_ext = [element.split('.')[1] for element in interface_args.filename]
if len(interface_args.filename) == 0:
    interface.print_help()
    print('='*42)
    sys.exit('Solution: Please consider inputting files' + '\n' + '='*42)

elif len(file_list_ext) != len(set(file_list_ext)):
    print('=' * 42)
    sys.exit('Duplicate file formats detected, only insert one of each file types' + '\n' + '=' * 42)

elif set(file_list_ext) != {'gff', 'fasta', 'vcf'}:
    print('=' * 42)
    missing = [ext for ext in ['fasta', 'gff', 'vcf'] if ext not in file_list_ext]
    for missing in missing:
        print(f'Missing {missing} file')
    sys.exit(f'Exit due to missing above files' + '\n' + '=' * 42)

else:
    for file in interface_args.filename:
        if file.endswith(('.vcf', 'vcf.gz')):
            vcf_data = file
        elif file.endswith('.gff'):
            gff_data = file
        elif file.endswith('.fasta'):
            fasta_data = file
        else:
            interface.print_help()
            print('=' * 42)
            sys.exit(f'{file} is not an acceptable file format' + '\n' + '=' * 42)

        file_list_str += file.split("/")[-1] + ', '

logging.info(f'Files inputted: {file_list_str[:-2]}')

# Directory
path = os.getcwd().split('/', maxsplit=(os.getcwd().count('/') - 1))[-1]

warnings.simplefilter('ignore', Bio.BiopythonWarning)

# VCF parser
reader = vcf.VCFReader(filename=vcf_data)

# GFF DB
try:
    gff_db = gffutils.create_db(gff_data, dbfn='PlasmoDB-54.db')
    logging.info(f'Creating PlasmoDB-54.db from {gff_data}')
except sqlite3.OperationalError:  # Passing the step if the DB exists
    gff_db = gffutils.FeatureDB('PlasmoDB-54.db', keep_order=True)
    logging.info(f'Accessing PlasmoDB-54.db')


# This function takes both the reference and alt codon of each SNP and check if the mutation is silent or results in
# an amino acid change, returns a binary output for the counter variable and a text output for the results tsv file
def mutation_qualifier(ref, alt):
    try:
        if Seq.translate(ref) == Seq.translate(alt):
            mutation_type = 'Synonymous'
            return 1, mutation_type
        else:
            mutation_type = 'Non-synonymous'
            return 0, mutation_type
    except KeyError:
        logging.info('Error at this line: ', line.CHROM, line.POS, ref, alt)


# Counters
below_quality = 0
above_quality = 0
synonymous_count = 0
non_syno_count = 0

# DF initialiser
tsv_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'Type',
                               'Transcript', 'Protein Location',
                               'Ref AA', 'Alt AA'])

# VCF parser
logging.info('Parsing VCF file')
print('Starting iteration through VCF')
for line in reader:
    if line.QUAL > 20:
        above_quality += 1
        # print(line.CHROM, line.POS)
        # Filtering for coding region only
        for region in gff_db.region(seqid=line.CHROM, featuretype='CDS', start=line.start, end=line.end):
            seq = ''  # Initialise seq string
            for parent in gff_db.parents(id=region.id, featuretype='mRNA'):
                mutation_position = 0  # Initialise mutation coordinates
                if parent.strand == '+':
                    # Iterating through all children to get full sequence and mutation coords
                    for child in gff_db.children(parent.id, featuretype='CDS', reverse=False, order_by='start'):
                        if child.attributes['Parent'][0] == parent.id:  # Checking if GFF file is correct
                            seq += child.sequence(fasta=fasta_data, use_strand=True)  # strings every cds seq together

                        else:  # Gff file error handler
                            logging.info(f"For {line.CHROM} {line.POS}, the childs' parent ID and the"
                                         f"parent ID are not the same. {child.attributes['Parent'][0]} is not "
                                         f"{parent.id}, potential error with GFF file, please inspect. Script will"
                                         f" terminate after this message")
                            logging.shutdown()
                            continue

                    for child in gff_db.children(parent.id, featuretype='CDS', reverse=False, order_by='start'):
                        if child.attributes['Parent'][0] == parent.id:
                            if child.start < line.POS < child.end:
                                # Here we select children CDS from mRNA parent where the mutation lies between its range
                                mutation_position += int(line.start) - child.start + 1  # + 1 for python counting
                                break
                            elif not child.start < line.POS < child.end:
                                # This corrects the mutation coords if parent contains multiple CDSs
                                mutation_position += child.end - child.start + 1  # + 1 to account for python counting
                            else:
                                print(f'Issue with this CDS: {child.id}')



                elif parent.strand == '-':  # The same is repeated for negative strand
                    for child in gff_db.children(parent.id, featuretype='CDS', reverse=True, order_by='start'):
                        if child.attributes['Parent'][0] == parent.id:  # Checking if GFF file is correct
                            seq += child.sequence(fasta=fasta_data, use_strand=True)  # strings every cds seq together
                        else:
                            logging.info(f"For {line.CHROM} {line.POS}, the childs' parent ID and the"
                                         f"parent ID are not the same. {child.attributes['Parent'][0]} is not "
                                         f"{parent.id}, potential error with GFF file, please inspect. Script will"
                                         f" terminate after this message")
                            logging.shutdown()
                            sys.exit(1)

                    for child in gff_db.children(parent.id, featuretype='CDS', reverse=True, order_by='start'):
                        if child.start < line.POS < child.end:
                            mutation_position += child.end - int(line.end)
                            break
                        else:
                            mutation_position += child.end - child.start


            for parent in gff_db.parents(id=region.id, featuretype='pseudogenic_transcript'):
                    seq = ''
                    logging.info(f'CDS: {region.id} is a pseudogenic transcript')
                    continue

            # Incorrect seq length handler
            #if len(seq) % 3 != 0:
               # logging.info(f'Sequence length for {region.id} not multiple of three, {len(seq) % 3} extra')

            # Generating protein sequence, calculating mutation position, empty seq handler
            try:
                prot = Seq(seq).translate()
                print(prot)
                protein_pos = math.ceil(mutation_position / 3)
            except NameError:
                if seq == '' and parent.featuretype != 'pseudogenic_transcript':
                    logging.info('Seq empty due to error with GFF file')
                elif seq == '' and parent.featuretype == 'pseudogenic_transcript':
                    continue
                break

            # Complement alt base pair based on strand
            if region.strand == '-':
                alt_bp = Seq(str(line.ALT[0])).reverse_complement()
                protein_pos += 1
            elif region.strand == '+':
                alt_bp = str(line.ALT[0])
            else:
                logging.info(f'No strand sense data for {region.id}')

            if prot.startswith('M') and prot[-1] == '*' and prot.count('*') == 1:  # Correct Protein seq test
                if mutation_position % 3 == 0:  # mutation at beginning of the codon
                    ref_codon = seq[mutation_position: mutation_position + 3]
                    alt_codon = alt_bp + ref_codon[1:]
                    if mutation_qualifier(ref_codon, alt_codon)[0] == 1:
                        synonymous_count += 1
                    else:
                        non_syno_count += 1

                elif mutation_position % 3 == 1:  # mutation in the middle of the codon
                    ref_codon = seq[mutation_position - 1: mutation_position + 2]
                    alt_codon = ref_codon[0] + alt_bp + ref_codon[2]
                    if mutation_qualifier(ref_codon, alt_codon)[0] == 1:
                        synonymous_count += 1
                    else:
                        non_syno_count += 1

                elif mutation_position % 3 == 2:  # mutation at the end of the codon
                    ref_codon = seq[mutation_position - 2: mutation_position + 1]
                    alt_codon = ref_codon[:2] + alt_bp
                    if mutation_qualifier(ref_codon, alt_codon)[0] == 1:
                        synonymous_count += 1
                    else:
                        non_syno_count += 1

            else:  # Incorrect protein seq handler
                logging.info(f'CDS: {region.id} for SNP {line.CHROM} {line.POS} '
                             f'has a pseudo-genetic protein sequence and was skipped')
                if not prot.endswith('*'):
                    logging.info('Protein seq did not have * at end')
                if not prot.startswith('M'):
                    logging.info('No Meth at start')
                continue

            # Translating codons and setting alt to NA if synonymous
            ref_aa = Seq.translate(Seq(ref_codon))
            alt_aa = Seq.translate(alt_codon)
            if ref_aa == alt_aa:
                alt_aa = 'NA'

            try:
                tsv_df.loc[len(tsv_df)] = [line.CHROM, line.POS, line.REF[0], line.ALT[0],
                                           mutation_qualifier(ref_codon, alt_codon)[1],
                                           region.id, protein_pos, ref_aa, alt_aa]
            except Exception:
                print('Error when appending results to pandas DF for:', line.CHROM, line.POS)
                continue

        # Append row to DF for non-coding mutations, IF-ELSE block handles duplicated entries error I was getting
        # saying that each VCF row had both non-coding and coding mutations
        if tsv_df.empty:
            tsv_df.loc[len(tsv_df)] = [line.CHROM, line.POS, line.REF[0], line.ALT[0],
                                       'Non-coding', 'NA', 'NA', 'NA', 'NA']
        elif not tsv_df.empty:
            if tsv_df.iloc[-1, 1] != line.POS and tsv_df.iloc[-1, 2] != line.REF and tsv_df.iloc[-1, 3] != line.ALT:
                tsv_df.loc[len(tsv_df)] = [line.CHROM, line.POS, line.REF[0], line.ALT[0],
                                           'Non-coding', 'NA', 'NA', 'NA', 'NA']

    elif line.QUAL <= 20:
        below_quality += 1
        continue
    else:
        logging.info(f'Error: {line.CHROM} {line.POS}, the quality: {line.QUAL} was erroneous')
        continue

# Final logs and TSV output
logging.info(f'{below_quality} VCF SNPs had read qualities below threshold of 20')
tsv_df.to_csv(sep='\t', path_or_buf='2117832_results.tsv', index=False)
logging.info(f'Final Results Output generated and saved as 2117832_results.tsv at {path}')
print('TSV outputted')

non_coding = above_quality - (synonymous_count + non_syno_count)
total = non_coding + synonymous_count + non_syno_count

# proportions
non_coding = (non_coding/total) * 100
synonymous_count = (synonymous_count/total) * 100
non_syno_count = (non_syno_count/total) * 100

# Bar plot generation
sns.barplot(
    x=['Non-coding', 'Synonymous', 'Non-synonymous'],
    y=[non_coding, synonymous_count, non_syno_count]
).set(
    title='Effect of Variant mutations',
    ylabel='Counts (%)'
)
logging.info(f'Bar plot created and outputted as 2117832_barplot.png saved at {path}')
plot.savefig('2117832_barplot.png')