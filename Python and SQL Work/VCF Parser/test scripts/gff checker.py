# IMPORTS
import vcf
import gffutils
from gffutils import biopython_integration
import csv
import os
import sqlite3
import seaborn as sns
import matplotlib.pyplot as plot
import pandas as pd
from Bio.Seq import Seq
import math

# DATA
toy_data = 'Toy Data/testData.vcf'
real_data = 'Assessment Data/assessmentData.vcf.gz'
gff_data = 'Genome files/PlasmoDB-54_Pfalciparum3D7.gff'
fasta_data = 'Genome files/PlasmoDB-54_Pfalciparum3D7_Genome.fasta'

# VCF parser
# reading toy data
reader = vcf.VCFReader(filename=toy_data)

# GFF DB
try:
    gff_db = gffutils.create_db(gff_data, dbfn='PlasmoDB-54.db')

# Passing the step if the DB exists
except sqlite3.OperationalError:
    gff_db = gffutils.FeatureDB('PlasmoDB-54.db', keep_order=True)


def mutation_qualifier(ref, alt):
    try:
        if Seq.translate(ref) == Seq.translate(alt):

            mutation_type = 'Non-synonymous'
            return 1, mutation_type
        else:

            mutation_type = 'Synonymous'
            return 0, mutation_type
    except KeyError:
        print('Error at this line: ', line.CHROM, line.POS, ref, alt)


row_counter = 0
coding_variant = 0
synonymous_count = 0
non_syno_count = 0

# DF initialiser
tsv_df = pd.DataFrame(columns=['CHROM', 'POS', 'REF', 'ALT', 'Type',
                               'Transcript', 'Protein Location',
                               'Ref AA', 'Alt AA'])
for line in reader:
    if line.QUAL > 20:
        print('#' * 30)
        print(line.CHROM, line.POS)
        print('#' * 30)
        row_counter += 1
        for region in gff_db.region(seqid=line.CHROM, featuretype='CDS', start=line.start, end=line.end):
            for parent in gff_db.parents(id=region.id, featuretype='mRNA'):
                mutation_position = 0  # Initialise mutation coordinates
                seq = ''  # Initialise seq string
                if parent.strand == '+':
                    # Filtering for vcf_rows region only
                    for child in gff_db.children(parent.id, featuretype='CDS', reverse=False, order_by='start'):
                        if child.attributes['Parent'][0] == parent.id:  # Logic test to
                            seq += child.sequence(fasta=fasta_data, use_strand=True)  # strings every cds seq together
                    print(seq)

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

                # The same is repeated for negative strand
                elif parent.strand == '-':
                    for child in gff_db.children(parent.id, featuretype='CDS', reverse=True, order_by='start'):
                        seq += child.sequence(fasta=fasta_data, use_strand=True)

                    for child in gff_db.children(parent.id, featuretype='CDS', reverse=False, order_by='start'):
                        if child.attributes['Parent'][0] == parent.id:
                            if child.start < line.POS < child.end:
                                mutation_position += child.end - int(line.end)
                                break
                            else:
                                mutation_position += child.end - child.start

            coding_variant += 1
            #
            prot = Seq(seq).translate()
            protein_pos = math.ceil(mutation_position / 3)

            # Alt base pair based on strand
            if region.strand == '-':
                alt_bp = Seq(str(line.ALT[0])).reverse_complement()
                protein_pos += 1

            elif region.strand == '+':
                alt_bp = str(line.ALT[0])

            else:
                print(f'No strand sense for {region.id}')

            print(protein_pos)

            if len(seq) % 3 != 0:
                print('Sequence length not multiple of three,', len(seq) % 3, 'extra')

            if prot.startswith('M') and prot[-1] == '*' and prot.count('*') == 1:  # Protein seq verifier
                if mutation_position % 3 == 0:
                    print('mutation at beginning of the codon')
                    ref_codon = seq[mutation_position: mutation_position + 3]
                    print(f'REF codon: {ref_codon}, {Seq.translate(ref_codon)}')
                    print(f'REF: {line.REF[0]} ALT: {line.ALT[0]}')
                    alt_codon = alt_bp + ref_codon[1:]
                    print(f'alt codon: {alt_codon}, {Seq.translate(alt_codon)}')

                    if mutation_qualifier(ref_codon, alt_codon) == 1:
                        synonymous_count += 1
                    else:
                        non_syno_count += 1

                elif mutation_position % 3 == 1:
                    print('mutation in the middle of the codon')
                    ref_codon = seq[mutation_position - 1: mutation_position + 2]
                    print(f'REF codon: {ref_codon}, {Seq.translate(ref_codon)}')
                    print(f'REF: {line.REF[0]} ALT: {line.ALT[0]}')
                    alt_codon = ref_codon[0] + alt_bp + ref_codon[2]
                    print(f'alt codon: {alt_codon}, {Seq.translate(alt_codon)}')

                    if mutation_qualifier(ref_codon, alt_codon) == 1:
                        synonymous_count += 1
                    else:
                        non_syno_count += 1

                elif mutation_position % 3 == 2:
                    print('mutation at the end of the codon')
                    ref_codon = seq[mutation_position - 2: mutation_position + 1]
                    print(f'REF codon: {ref_codon}, {ref_codon.translate()}')
                    print(f'REF: {line.REF[0]} ALT: {line.ALT[0]}')
                    print('region mutation bp', seq[mutation_position + 1])
                    alt_codon = ref_codon[:2] + alt_bp
                    print(f'alt codon: {alt_codon}, {Seq.translate(alt_codon)}')

                    if mutation_qualifier(ref_codon, alt_codon)[0] == 1:
                        synonymous_count += 1

                    else:
                        non_syno_count += 1


            else:
                pass
                # Log error in logger

    if row_counter > 0:
        continue


# Finds all ACT codons in gene seq
# act = [index for index in range(len(seq)) if seq.startswith('ACT', index)]
