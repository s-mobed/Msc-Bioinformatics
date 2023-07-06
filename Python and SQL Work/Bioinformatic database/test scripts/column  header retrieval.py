# Imports
from contextlib import closing
import sqlite3
import csv
import glob
datafiles = ['HMP_metabolome_abundance.tsv', 'HMP_metabolome_annotation.csv',
             'HMP_proteome_abundance.tsv', 'HMP_transcriptome_abundance.tsv',
             'Subject.csv']
for file in datafiles:
    with open(file) as open_file:
        print(file.rsplit(sep='.')[0])
        for line in open_file:
            header_list = line.replace(',', ' varchar(255),').replace('\t', ' varchar(255),')\
                .replace('\n', ' varchar(255)').replace('\ufeff', '')
            print(header_list)
            break
