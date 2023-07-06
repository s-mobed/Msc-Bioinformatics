#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=24:00:00
#PBS -N caller
#PBS -d /export/home/biostuds/2117832m/BIOL5299/
#PBS -m abe
#PBS -M 2117832m@student.gla.ac.uk
#PBS -q bioinf-stud




ref=/export/home/biostuds/2117832m/BIOL5299/raw_data/Reference/TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa

mkdir -p VCF

bamaddrg -b bowtie_results/WT.bam -s WT -b bowtie_results/AmpB.bam -s AmpB | freebayes -f $ref -p 2 --stdin > VCF/var.vcf

vcffilter -f "QUAL < 20 & TYPE = snp" VCF/var.vcf > VCF/alt_filtered.vcf

picard CreateSequenceDictionary R=$ref O=Lmex.dict


