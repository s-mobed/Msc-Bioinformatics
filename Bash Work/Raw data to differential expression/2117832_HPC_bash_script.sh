#PBS -l nodes=1:ppn=4:centos6,cput=24:00:00,walltime=48:00:00 
#PBS -N assessment_1 
#PBS -d /export/home/biostuds/2117832m/BIOL5177/part1/ 
#PBS -m abe 
#PBS -M 2117832m@student.gla.ac.uk 
#PBS -q bioinf-stud 

# 
# RESOURCE FILES 

adapter=/export/projects/polyomics/biostuds/data/illumina_adapter.fa # path to illumina adapter file 
hs2index=/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2 # path to reference genome hisat2 indexes 
gtf=/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf  # path to GTF file 
data=/export/projects/polyomics/buzz/biostuds # path to local data directory that you have created 

# 
# MAKE FEW SUBDIRS UNLESS THEY EXIST 

mkdir -p hisat_results # make above hisat2 directory, -p checks if such dir already exists 

mkdir -p stringtie_results # make above stringtie directory 

gtflist='list.gtf.txt' # filename for final GTF list 
rm -f ${gtflist} # remove if exists 


#RUNNING a single LOOP for all the work 

for sample in s1.c2 s2.c2 s3.c2 s4.c2 s5.c2 s6.c2 s7.c2 s8.c2 s9.c2 s10.c2 s11.c2 s12.c2
do 

# Setting temporary variables

	fastq="$data/${sample}.fq" # path to raw fastq file 
	trim1="${sample}.t1.fq" # path to adapter-trimmed fastq file 
	trim2="${sample}.t2.fq" # path to quality-trimmed fastq file 

# Adapter and Quality Trimming

	scythe $fastq -o $trim1 -a $adapter -q sanger

	sickle se -f $trim1 -o $trim2 -q 10 -t sanger -l 49

# Alignment of single end reads

	hisat2 -p 4 --rna-strandness 'R'  --phred33 -x $hs2index -U $trim2 -S hisat_results/${sample}.sam

# Sam to bam conversion through direct pipe to avoid creating extra variable, and quicker processing

	samtools view -@ 4 -Sbu -o - hisat_results/${sample}.sam | samtools sort -@ 4 - hisat_results/${sample}.sort 

	
# removing unnecessary files

	rm hisat_results/${sample}.sam # removing unnecessary files 
	rm ${trim1} ${trim2} # removing unnecessary files 

# creating sample specific sub directories

	str_smp_dir="stringtie_results/${sample}-tie" # path to sample-specific subdirectory for stringtie results 

	mkdir -p $str_smp_dir # make the above directory 


# transcript assembly and quantification

	stringtie -e -B -p 4 -G $gtf -o $str_smp_dir/${sample}.gtf hisat_results/${sample}.sort.bam

	gtfline="${sample} $str_smp_dir/${sample}.gtf" # line containing sample and path to GTF

	echo ${gtfline} >> ${gtflist} # adding the above line to a file

done

# Converting sample-specific GTFs to a single gene-count matrix
python2.7 /export/projects/polyomics/App/prepDE.py -i ${gtflist} -g gene_count.csv
