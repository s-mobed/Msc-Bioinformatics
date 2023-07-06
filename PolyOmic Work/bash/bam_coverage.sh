#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=24:00:00
#PBS -N caller
#PBS -d /export/home/biostuds/2117832m/BIOL5299/
#PBS -m abe
#PBS -M 2117832m@student.gla.ac.uk
#PBS -q bioinf-stud




ref=/export/home/biostuds/2117832m/BIOL5299/raw_data/Reference/TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa

for sample in AmpB WT
do

	bamCoverage -b bowtie_results/${sample}.bam -o coverage/${sample}.bed

done
