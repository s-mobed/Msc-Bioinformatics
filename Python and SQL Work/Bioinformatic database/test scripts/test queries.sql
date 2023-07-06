-- 1;
SELECT SubjectID, Age FROM Subject WHERE Age  > 70;

-- 2;
SELECT SubjectID, BMI, Sex FROM Subject WHERE Sex = 'F' AND BMI  BETWEEN 18.5 AND 24.9 ORDER BY BMI DESC;

-- 3;
SELECT substr(SampleID,instr(SampleID,'-')+1) AS 'ZNQOVZV visits' FROM HMP_metabolome_abundance WHERE SampleID  LIKE 'ZNQOVZV%';

-- 4;
SELECT DISTINCT substr(HMP_metabolome_abundance.SampleID, 1, instr(SampleID, '-' ) -1 )  AS Sample_names, Subject.IR_IS_classification 
FROM HMP_metabolome_abundance, Subject
WHERE Subject.SubjectID = Sample_names AND IR_IS_classification = 'IR';

-- 5;
SELECT DISTINCT KEGG FROM HMP_metabolome_annotation WHERE
PeakID LIKE '%nHILIC_121.0505_3.5%' OR PeakID LIKE '%nHILIC_130.0872_6.3%' OR
PeakID LIKE '%nHILIC_133.0506_2.3%' OR PeakID LIKE '%nHILIC_133.0506_4.4%';

--6;
SELECT AVG(Age) AS Mean_Age, min(Age) AS Lowest_Age, max(Age) AS Highest_Age FROM Subject ;

--7;
SELECT Pathway, count(Pathway) AS Pathway_count FROM HMP_metabolome_annotation GROUP BY Pathway HAVING Pathway_count > 10 ORDER BY Pathway_count DESC;

--8;
SELECT max(A1BG) FROM HMP_transcriptome_abundance WHERE SampleID LIKE 'ZOZOW1T%';

--9;
SELECT BMI, AGE from Subject WHERE AGE AND BMI NOT NULL


