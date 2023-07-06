from contextlib import closing
import sqlite3
import pandas as pd
import seaborn as sns



def query_db(db_path, query_number):
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
    query_number = int(query_number)
    if query_number not in range(1, 10):
        return print('Error: Choose query number between 1 and 9 only')

    # Query closing object
    with closing(sqlite3.connect(db_path)) as connection:
        with closing(connection.cursor()) as cur:
            sql = sql_options[query_number]
            print(sql)
            results_list = []
            for results in cur.execute(sql):
                print('')
                for x in range(len(results)):
                    print(f"{results[x]}", end='\t')

                if query_number == 9:
                    results_list.append(results)

            # Creation of BMI-AGE dataframe
            plot_df = pd.DataFrame(results_list)
            plot_df = plot_df.rename(columns={0: 'BMI', 1: 'AGE'})  # Renaming columns
            plot_df = plot_df.astype(float)  # converting strings to floats

            # AGE-BMI plot generation
            plot = plot_df.plot(
                x='AGE',
                y='BMI',
                kind='scatter',
                title='Age-BMI relation from multi-omics study of aging \n n= 103, DOI: 10.1038/s41591-019-0719-5',

            )
            # Saving Figure & adding trend line
            sns.regplot(plot_df, x='AGE', y='BMI').get_figure().savefig(fname='age_bmi_scatterplot.png')
            print('')

db_path = '/Users/Shean/Desktop/MSc Bioinformatics/Programming and Databases/' \
          'Python labs/BIOL4292/Assignment 2/Assignment.db'

query_db(db_path, 9)
