# Imports
import csv
from contextlib import closing
import sqlite3
from collections import defaultdict
import re
from itertools import islice
from annotation_parser import annotation_parser

db_path = '/Users/Shean/Desktop/MSc Bioinformatics/Programming and Databases/' \
          'Python labs/BIOL4292/Assignment 2/Assignment.db'

datafiles = ['HMP_metabolome_abundance.tsv', 'HMP_metabolome_annotation.csv',
             'HMP_proteome_abundance.tsv', 'HMP_transcriptome_abundance.tsv',
             'Subject.csv']

for file in datafiles:
    with open(file) as open_file:
        file_name = file.rsplit(sep='.')[0]
        print(file.rsplit(sep='.')[0])
        # File type handler if/else statement
        if '.tsv' in file:
            reader = csv.reader(open_file, delimiter='\t')
        else:
            reader = csv.reader(open_file, delimiter=',')

        # Header and qmark question mark inserter
        header = [cell.replace('\ufeff', '') for cell in next(reader)]  # Stores and skips first row
        inserter_size = '?, ' * (len(header) - 1) + '?'

        # Separate data loading for annotation file
        if 'annotation' in file_name:
            annotation_parser(reader, db_path, file_name, header)

        else:
            for line in reader:
                for cell in line:
                    # None type standardised for all tables
                    line = [None if cell in ('', 'NA', 'Na', 'Unknown', 'unknown') else cell for cell in line]

                    # Unwanted keyword stripper
                    if cell.startswith('\ufeff'):
                        line = [cell.replace('\ufeff', '') for cell in line]

                # Context manager for data loading SQL statements
                with closing(sqlite3.connect(db_path, )) as connection:
                    with closing(connection.cursor()) as cur:
                        '''
                        Statement that deals with files with columns numbers > 999, got an operational error
                        where max number of variables were used up in the SQL qmark style for both abundance
                        files with too large a header. Didn't want to change fixed value as I would needed to 
                        recompile sqlite3. 
                        '''
                        # For small headers
                        if len(header) < 999:  # MAX_VARIABLE_NUMBER
                            sql = 'INSERT INTO ' + f'{file_name} ' + f'VALUES({inserter_size});'
                            cur.execute(sql, line)

                        # For large headers, doesn't insert NULL values properly however abundance files do not have
                        # null values present, but doesn't work well if other large files have null values present.
                        else:
                            sql = 'INSERT INTO ' + f'{file_name}{tuple(header)} ' + f'VALUES {tuple(line)};'
                            cur.execute(sql)

                    connection.commit()
