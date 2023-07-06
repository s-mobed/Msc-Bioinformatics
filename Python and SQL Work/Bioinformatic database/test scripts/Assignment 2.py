# Imports
from contextlib import closing
import sqlite3
import csv
from annotation_parser import annotation_parser
from collections import defaultdict
import re
import pandas as pd
import seaborn as sns
import argparse
from sys import argv

'''
This program allows for information storage and retrieval into a SQL database.This program comprises of three functions,
database creation, data loading, and information retrieval based on a set of nine specific queries.
The data in question stems from a multi-omics study comprised five data files.
'''

# Generating command line interface
db_parser = argparse.ArgumentParser(description='Data file parser, database creation and database loader functions')

db_parser.add_argument('--createdb', required=False, action='store_true',
                       help='Creates a db table from each file and uses first row as table header')

db_parser.add_argument('--loaddb', required=False, action='store_true',
                       help='Loads data into a db table that is created beforehand with --createdb')

db_parser.add_argument('--querydb', type=int, required=False,
                       help='Creates a db table from each file and uses first row as table header')

db_parser.add_argument('database_files', type=str, nargs='*', help='Input file(s) desired')

db_args = db_parser.parse_args()


class db_load_and_fetch:
    def __init__(self):
        self.file_list = []

        self.query_number = db_args.querydb

        self.options_check = len(argv)  # This checks how many options have been chosen

    '''This handles file inputs from the commandline, raises error if file not csv/tsv'''

    def input_handler(self, input):
        try:
            if type(input) is not list:
                raise TypeError
        except TypeError:
            print('Input should be list')

        for file in input:
            try:
                if file.endswith('.csv') or file.endswith('.tsv'):
                    self.file_list.append(file)
                else:
                    raise NameError

            except NameError:
                print(f'{file} is incorrect file type consider using .csv/tsv files')

    '''This creates tables using the information from the first row of the file as table headers'''

    def create_db(self, go):
        if go is True:
            # For loop that creates custom SQL create table statement per data file
            for file in self.file_list:
                with open(file) as open_file:
                    file_name = file.rsplit(sep='.')[0]  # Filename extraction

                    for line in open_file:

                        # Setting header for each table
                        header = '"' + line.replace(',', '" varchar(255),"').replace('\t', '" varchar(255),"') \
                            .replace('\n', '" varchar(255),"').replace('\ufeff', '')
                        header = header[:len(header) - 2]

                        # Primary key designation
                        # Setting Metabolite column as key for annotations
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

                        break
                    # Cursor object and SQL create statement execution
                    with closing(sqlite3.connect('Assignment.db')) as connection:
                        with closing(connection.cursor()) as cur:
                            try:
                                sql = 'CREATE TABLE ' + f'{file_name}( {header}' + f',PRIMARY KEY({unique_key}));'
                                cur.executescript(sql)


                            except sqlite3.OperationalError:
                                print(f'It appears that a table named {file_name} is already present\n'
                                      f'in the database. The programme uses the filename as the table name\n'
                                      f'in that case change the filename and try again if this is a new dataset\n')
                                continue


    def loader(self, go):
        if go is True:

            for file in self.file_list:
                with open(file) as open_file:
                    file_name = file.rsplit(sep='.')[0]

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
                        annotation_parser(reader, 'Assignment.db', file_name, header)

                    else:
                        for line in reader:
                            for cell in line:
                                # None type standardised for all tables
                                line = [None if cell in ('', 'NA', 'Na', 'Unknown', 'unknown')
                                        else cell for cell in line]

                                # Unwanted keyword stripper
                                if cell.startswith('\ufeff'):
                                    line = [cell.replace('\ufeff', '') for cell in line]

                            # Context manager for data loading SQL statements
                            with closing(sqlite3.connect('Assignment.db')) as connection:
                                with closing(connection.cursor()) as cur:
                                    try:
                                        '''
                                        Statement that deals with files with columns numbers > 999, got an operational
                                        error where max number of variables were used up in the SQL qmark style for
                                        both abundance files with too large a header. Didn't want to change fixed
                                        value as I would needed to recompile sqlite3. 
                                        '''
                                        # For small headers
                                        if len(header) < 999:  # MAX_VARIABLE_NUMBER
                                            sql = 'INSERT INTO ' + f'{file_name} ' + f'VALUES({inserter_size});'
                                            cur.execute(sql, line)

                                            ''' For large headers, doesn't insert NULL values properly however abundance
                                            files do not have null values present, but doesn't work well if
                                            other large files have null values present.
                                            '''
                                        else:
                                            sql = 'INSERT INTO ' + f'{file_name}{tuple(header)} ' + f'VALUES {tuple(line)};'
                                            cur.execute(sql)
                                    except sqlite3.IntegrityError:
                                        print(f'The SQL statement: {sql, line} returned an integrity error')
                                        continue
                                    except sqlite3.OperationalError:
                                        print('sqlite3.OperationalError occurred check if the table of the data file '
                                              'has been created first with the --createdb option')
                                        break

                                connection.commit()

    '''
    Created separate static method that parses and load the data from the annotation file for better readability
    '''

    @staticmethod
    def annotation_parser(reader, db_path, file_name, header):
        # Importing data from main load module
        db_path = db_path
        file_name = file_name
        # Peak id default dict list
        peak_id = defaultdict(list)

        for line in reader:
            # parses multiple compound entries
            if '|' in line[1] and not re.search(pattern=r"\((\d+)\)", string=line[1]):

                # Splitting multiple compounds
                split_name = line[1].split(sep='|')
                split_KEGG = line[2].split(sep='|')
                split_HMDB = line[3].split(sep='|')

                # Filling peak_id dict with PeakIDs as values and metabolites as keys
                for name in split_name:
                    name = name.strip()  # strips any whitespace around names after split
                    if name in peak_id:
                        peak_id[name].append(line[0])
                    else:
                        peak_id[name] = [line[0]]

                # Filling Empty KEGG and HMDB
                if line[2] == '':
                    split_KEGG = [None, None]

                if line[3] == '':
                    split_HMDB = [None, None]

                # creating matrix of split lines
                split_lines = [
                    (line[0], split_name[0].strip(), split_KEGG[0], split_HMDB[0], line[4], line[5]),
                    (line[0], split_name[1].strip(), split_KEGG[1], split_HMDB[1], line[4], line[5])
                ]

                # Annotations table inserter context manager
                with closing(sqlite3.connect(db_path)) as connection:
                    with closing(connection.cursor()) as cur:
                        try:
                            for item in split_lines:
                                sql = 'INSERT INTO ' + f'{file_name}{tuple(header)} ' + f'VALUES(?, ?, ?, ?, ?, ?);'
                                cur.execute(sql, item)

                        # Duplicate of Pro-Cys caused integrity error, this error handle creates an
                        # UPDATE command that goes over the duplicates and fills in the missing values
                        # and appends the duplicates peakID into the peakID column
                        except sqlite3.IntegrityError:
                            for item in split_lines:
                                sql = f'UPDATE {file_name} SET KEGG = ?, HMDB = ?, PeakID = ? WHERE Metabolite = ?;'
                                params = (item[2], item[3], str(tuple(peak_id[item[1]])), item[1])
                                cur.execute(sql, params)

                    connection.commit()

            # parses multiple peaks entries
            elif '|' not in line[1] and re.search(pattern=r"\((\d+)\)", string=line[1]):
                # Filling empty cells with None
                line = [None if cell in ('', 'NA', 'Na', 'Unknown', 'unknown') else cell for cell in line]

                # Strips end numbers
                numberless_metabolite = line[1][:-3]

                # appending new PeakIds to dict
                if numberless_metabolite in peak_id:
                    peak_id[numberless_metabolite].append(line[0])
                else:
                    peak_id[numberless_metabolite] = [line[0]]

                # Modified line to include multiple peak ids
                new_line = (
                    str(tuple(peak_id[numberless_metabolite])),
                    numberless_metabolite,
                    line[2], line[3], line[4], line[5]
                )
                # Annotations table inserter context manager
                with closing(sqlite3.connect(db_path)) as connection:
                    with closing(connection.cursor()) as cur:
                        try:
                            sql = 'INSERT INTO ' + f'{file_name}{tuple(header)} ' + f'VALUES(?, ?, ?, ?, ?, ?);'
                            cur.execute(sql, new_line)

                        except sqlite3.IntegrityError:

                            sql = f'UPDATE {file_name} SET PeakID = ? WHERE Metabolite = ?;'
                            params = (new_line[0], new_line[1])
                            cur.execute(sql, params)

                        connection.commit()

            # parses bot multiple peaks and compounds entries
            elif '|' in line[1] and re.search(pattern=r"\((\d+)\)", string=line[1]):
                # Filling empty cells with None
                line = [None if cell in ('', 'NA', 'Na', 'Unknown', 'unknown') else cell for cell in line]

                # Removing numbers from metabolite names
                numberless_metabolite = line[1][:-3]

                # Splitting multiple compounds
                split_name = numberless_metabolite.split(sep=' | ')

                # Filling peak_id dict with PeakIDs as values and metabolites as keys
                for name in split_name:
                    if name in peak_id:
                        peak_id[name].append(line[0])
                    else:
                        peak_id[name] = [line[0]]

                # creating new lines from split metabolites
                split_lines = ([str(tuple(peak_id[split_name[0]])), split_name[0],
                                line[2], line[3], line[4], line[5]],
                               [str(tuple(peak_id[split_name[1]])), split_name[1],
                                line[2], line[3], line[4], line[5]]
                               )

                # Annotations table inserter context manager
                with closing(sqlite3.connect(db_path)) as connection:
                    with closing(connection.cursor()) as cur:
                        try:
                            for item in split_lines:
                                sql = 'INSERT INTO ' + f'{file_name}{tuple(header)} ' + f'VALUES(?, ?, ?, ?, ?, ?);'
                                cur.execute(sql, item)

                        except sqlite3.IntegrityError:
                            for item in split_lines:
                                sql = f'UPDATE {file_name} SET PeakID = ? WHERE Metabolite = ?;'
                                params = (item[0], item[1])
                                cur.execute(sql, params)

                        connection.commit()

            # Simple entries
            else:
                # Filling empty cells with None
                line = [None if cell in ('', 'NA', 'Na', 'Unknown', 'unknown') else cell for cell in line]

                # Annotations table inserter context manager
                with closing(sqlite3.connect(db_path)) as connection:
                    with closing(connection.cursor()) as cur:
                        try:
                            sql = 'INSERT INTO ' + f'{file_name}{tuple(header)} ' + f'VALUES(?, ?, ?, ?, ?, ?);'
                            cur.execute(sql, line)
                        except sqlite3.IntegrityError:
                            print('Statement not committed -->', line)
                            continue
                    connection.commit()


    def query_db(self, db_path, query_number):
        # Dictionary of SQL commands
        sql_options = {
            1: 'SELECT SubjectID, Age FROM Subject WHERE Age  > 70;',

            2: "SELECT SubjectID, BMI, Sex FROM Subject WHERE Sex = 'F' AND BMI  BETWEEN 18.5 AND 24.9 ORDER BY BMI DESC;",

            3: "SELECT substr(SampleID,instr(SampleID,'-')+1) AS 'ZNQOVZV visits' FROM HMP_metabolome_abundance WHERE"
               " SampleID  LIKE 'ZNQOVZV%';",

            4: "SELECT DISTINCT substr(HMP_metabolome_abundance.SampleID, 1, instr(SampleID, '-' ) -1 )  AS Sample_names,"
               " Subject.IR_IS_classification FROM HMP_metabolome_abundance, Subject"
               " WHERE Subject.SubjectID = Sample_names AND IR_IS_classification = 'IR';",

            5: "SELECT DISTINCT KEGG FROM HMP_metabolome_annotation WHERE PeakID LIKE '%nHILIC_121.0505_3.5%' "
               "OR PeakID LIKE '%nHILIC_130.0872_6.3%' OR"
               " PeakID LIKE '%nHILIC_133.0506_2.3%' OR PeakID LIKE '%nHILIC_133.0506_4.4%';",

            6: 'SELECT AVG(Age) AS Mean_Age, min(Age) AS Lowest_Age, max(Age) AS Highest_Age FROM Subject ;',

            7: 'SELECT Pathway, count(Pathway) AS Pathway_count FROM HMP_metabolome_annotation'
               ' GROUP BY Pathway HAVING Pathway_count > 10 ORDER BY Pathway_count DESC;',

            8: "SELECT max(A1BG) FROM HMP_transcriptome_abundance WHERE SampleID LIKE 'ZOZOW1T%';",

            9: 'SELECT BMI, AGE from Subject WHERE AGE AND BMI NOT NULL'
        }

        # Query number error handler
        if query_number is None:
            print('query is nonetype')
        elif type(query_number) is str:
            query_number = int(query_number)

        elif query_number not in range(1, 10):
            return print(f'Error with {query_number}: Choose query number between 1 and 9 only')

        # Query closing object
        with closing(sqlite3.connect(db_path)) as connection:
            with closing(connection.cursor()) as cur:
                try:
                    sql = sql_options[query_number]  # Retrieves sql command from dictionary
                    results_list = []
                    for results in cur.execute(sql):
                        print('')
                        for x in range(len(results)):
                            print(f"{results[x]}", end='\t')
                            results_list.append(results)

                    if len(results_list) == 0:
                        print(f'No results for statement:\n{sql}')

                    if query_number == 9:


                        # Creation of BMI-AGE dataframe
                        plot_df = pd.DataFrame(results_list)
                        plot_df = plot_df.rename(columns={0: 'BMI', 1: 'AGE'})  # Renaming columns
                        plot_df = plot_df.astype(float)  # converting strings to floats

                        # AGE-BMI plot generation
                        plot = plot_df.plot(
                            x='AGE',
                            y='BMI',
                            kind='scatter',
                            title='Age-BMI relation from multi-omics study of aging'
                                  ' \n n= 103, DOI: 10.1038/s41591-019-0719-5'

                        )
                        # Saving Figure & adding trend line
                        sns.regplot(data=plot_df, x='AGE', y='BMI').get_figure().\
                            savefig(fname='age_bmi_scatterplot.png')
                except sqlite3.OperationalError:
                    print(f'Operational Error occurred with: {sql}, perhaps the table has not been created'
                          f'or loaded with data yet.')
            print('')


    '''Initialise class'''


instance = db_load_and_fetch()
instance.input_handler(db_args.database_files)
if db_args.createdb:
    instance.create_db(db_args.createdb)

elif db_args.loaddb:
    instance.loader(db_args.loaddb)

elif instance.query_number is not None:
    instance.query_db('Assignment.db', db_args.querydb)

elif instance.options_check == 1:
    print('No options provided')
