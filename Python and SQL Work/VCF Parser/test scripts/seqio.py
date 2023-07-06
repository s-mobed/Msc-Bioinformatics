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
from Bio import SeqIO
import math

# DATA
toy_data = 'Toy Data/testData.vcf'
real_data = 'Assessment Data/assessmentData.vcf.gz'
gff_data = 'Genome files/PlasmoDB-54_Pfalciparum3D7.gff'
fasta_data = 'Genome files/PlasmoDB-54_Pfalciparum3D7_Genome.fasta'

# VCF parser
# reading toy data
reader = vcf.VCFReader(filename=toy_data)

for entry in SeqIO.parse(fasta_data, 'fasta'):
    if entry.id == 'Pf3D7_14_v3':
        seq = entry.seq[205558:208840]
        # seq = Seq(''.join([element for element in [str(seq)[codon:codon+3]for codon in range(0, len(seq), 3)]][::-1])).complement()
        print(seq)
        print(len(seq) % 3)
        prot = seq.translate()
        print(prot)

Seq(''.join([element for element in [str(seq)[codon:codon+3]for codon in range(0, len(seq), 3)]][::-1])).complement()