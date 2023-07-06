#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=24:00:00 
#PBS -N pe_pipeline
#PBS -d /export/home/biostuds/2117832m/BIOL5299/ 
#PBS -m abe 
#PBS -M 2117832m@student.gla.ac.uk 
#PBS -q bioinf-stud 

# 
# RESOURCE FILES 
wd=/export/home/biostuds/2117832m/BIOL5299/
index=/export/home/biostuds/2117832m/BIOL5299/bowtie_index/
data=/export/home/biostuds/2117832m/BIOL5299/raw_data/DNAseq # path to local data directory that you have created 

# 

#RUNNING a single LOOP for all the work 

for sample in LmexAmpB LmexWT
do 

# Setting var names for trimmed reads
	trim1="trimmed_results/${sample}_1_val_1.fq.gz" # path to adapter-trimmed fastq file 
	trim2="trimmed_results/${sample}_2_val_2.fq.gz" # path to quality-trimmed fastq file 
	
	echo $trim1
	echo $trim2
 
# Bowtie2 step
	#if [ -n "$(ls -A bowtie_results/)" ]; then
               #continue
	#else
	bowtie2 --phred64 -x $indexLmex -1 $trim1 -2 $trim2 -S bowtie_results/${sample}.sam
	#echo $bowtieCmd
	#fi

# Samtools step
	samtools view -@ 4 -Sbu -o - bowtie_results/${sample}.sam | samtools sort -@ 4 - bowtie_results/${sample}.sort

done 

