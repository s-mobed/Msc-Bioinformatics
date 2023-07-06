#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=24:00:00 
#PBS -N pe_pipeline
#PBS -d /export/home/biostuds/2117832m/BIOL5299/ 
#PBS -m abe 
#PBS -M 2117832m@student.gla.ac.uk 
#PBS -q bioinf-stud 

# 
# RESOURCE FILES 
wd=/export/home/biostuds/2117832m/BIOL5299/
index=/export/home/biostuds/2117832m/BIOL5299/bowtie_index/Lmex
data=/export/home/biostuds/2117832m/BIOL5299/raw_data/DNAseq # path to local data directory that you have created 

# 
# MAKE FEW SUBDIRS UNLESS THEY EXIST 

mkdir -p bowtie_results # make above bowtie directory, -p checks if such dir already exists 

mkdir -p bowtie_index # creating leishmania index 

mkdir -p fastqc_results # make above stringtie directory 

mkdir -p trimmed_results

# Bowtie build reference index

# bowtie2-build --threads 4 /export/home/biostuds/2117832m/BIOL5299/raw_data/Reference/*.fa bowtie_index/Lmex


#RUNNING a single LOOP for all the work 
echo 'before loop'
for sample in LmexAmpB LmexWT
do
        echo $sample

        # Setting var names for raw paired reads
        fastq1="$data/${sample}_1.fastq.gz" # path to raw fastq file
        fastq2="$data/${sample}_2.fastq.gz" # path to raw fastq file

        # Fastqc step
        if [ -f "fastqc_results/${sample}_1_fastqc.html" ] && [ -f "fastqc_results/${sample}_2_fastqc.html" ]; then
                echo 'fastqc already done'
        else
                fastqc -o fastqc_results/ -t 4 $fastq1 $fastq2
        fi

        # Setting var names for trimmed reads
        trim1="trimmed_results/${sample}_1_val_1.fq.gz" # path to adapter-trimmed fastq file
        trim2="trimmed_results/${sample}_2_val_2.fq.gz" # path to quality-trimmed fastq file

        # Trim Galore step
        if [ -f "trimmed_results/${sample}_1_val_1.fq.gz" ] && [ -f "trimmed_results/${sample}_2_val_2.fq.gz" ]; then
                echo 'trim already done'
        else
                trim_galore -q 20 --phred64 --illumina -o trimmed_results/ --paired $fastq1 $fastq2
        fi

        # Bowtie2 & samtools step
        if [ -f "bowtie_results/${sample}.sam" ]; then
                echo 'bowtie already done'
        else
		bowtie2 --phred64 -x $index -1 $trim1 -2 $trim2 -S $sample.sam
		
		samtools view -@ 4 -Sbu $sample.sam | samtools sort -@ 4 -o bowtie_results/${sample}.bam -
		
		samtools index bowtie_results/${sample}.bam
        fi

	bamCoverage -b $sample.bam -of 'bigwig' -o $sample.bw

	# Removing sam file
	rm ${sample}.sam
done

