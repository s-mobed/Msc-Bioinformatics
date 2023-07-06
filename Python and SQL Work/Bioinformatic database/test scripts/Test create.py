# Imports
from contextlib import closing
import sqlite3


datafiles = ['HMP_metabolome_abundance.tsv', 'HMP_metabolome_annotation.csv',
             'HMP_proteome_abundance.tsv', 'HMP_transcriptome_abundance.tsv',
             'Subject.csv']

# For loop that creates custom SQL create table statement per data file
for file in datafiles:
    with open(file) as open_file:
        # Filename extraction
        file_name = file.rsplit(sep='.')[0]
        # Checking output
        print(file.rsplit(sep='.')[0])
        for line in open_file:
            # Setting header for each table
            header = '"' + line.replace(',', '" varchar(255),"').replace('\t', '" varchar(255),"')\
                .replace('\n', '" varchar(255),"').replace('\ufeff', '')
            header = header[:len(header)-2]

            # Primary key designation
            # Setting Metabolite column as primary key
            if file_name == 'HMP_metabolome_annotation':
                if '\t' in line:
                    unique_key = line.split('\t')[1].replace('\ufeff', '')
                else:
                    unique_key = line.split(',')[1].replace('\ufeff', '')

            # Setting first column of the other tables as the key
            else:
                if '\t' in line:
                    unique_key = line.split('\t')[0].replace('\ufeff', '')
                else:
                    unique_key = line.split(',')[0].replace('\ufeff', '')

            # Cursor object and SQL create statement execution

            with closing(sqlite3.connect('Assignment.db')) as connection:
                with closing(connection.cursor()) as cur:
                    sql = 'CREATE TABLE ' + f'{file_name}( {header}' + f',PRIMARY KEY({unique_key}));'
                    print(sql)
                    cur.executescript(sql)

