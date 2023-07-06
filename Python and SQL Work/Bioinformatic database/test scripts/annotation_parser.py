# Imports
from contextlib import closing
import sqlite3
from collections import defaultdict
import re


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
            print(line)

            # Annotations table inserter context manager
            with closing(sqlite3.connect(db_path)) as connection:
                with closing(connection.cursor()) as cur:
                    try:
                        sql = 'INSERT INTO ' + f'{file_name}{tuple(header)} ' + f'VALUES(?, ?, ?, ?, ?, ?);'
                        cur.execute(sql, line)
                    except sqlite3.IntegrityError:
                        print('Faulty line -->', line)
                        continue
                connection.commit()
