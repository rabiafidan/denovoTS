#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID
#BSUB -J pipeline

#Data: High coverage 1KG trios
#https://www.internationalgenome.org/data-portal/data-collection/30x-grch38

#load r and python modules
ml load r-4.0.3-gcc-9.3.0-4l6eluj
ml load python-3.9.0-gcc-9.3.0-5t75egs

#create a directory for data manipulation
if [ ! -d data ]; then
    mkdir data;
fi

#create crams directory to store cram files
if [ ! -d crams ]; then
    mkdir crams;
fi

#download the reference fasta
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

#index the reference fasta
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

cd data

#fetch pedigree file
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

wd=$(pwd) # r/python working directory

#produce children.txt,father.txt,mother.txt,inds.txt,pops.txt and sup_pops.txt files
Rscript ../scripts/data.R $wd 

#download the table containing the cram file links from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
#You can also get the same list from the github repo by uncommenting the following line:
#wget https://raw.githubusercontent.com/rabiafidan/denovoTS/master/igsr_30x%20GRCh38.tsv?token=ANT5TMVRQM2OESEY3DFZAPLBIBS3O

#obtain cram file links of the all individuals of the trios
grep cram igsr_30x\ GRCh38.tsv | cut -f1 >links.txt
for i in $(cat inds.txt); do grep $i links.txt >> links2.txt; done
mv links2.txt links.txt

cd ..

#download the individual cram files
bsub -J "dwnld[1-1793]%40" < scripts/download.sh

#create the configuration file for Snakefile
python scripts/configfile.py $wd 

#create a bed file to use for regional filtering of vcf files
#Download UCSC Dec2013 hg38 repeatmasker file as repeatmask.txt. Selected fields (GenoName GenoStart GenoEnd repClass)
grep 'Low_complexity' repeatmask.txt  | grep $'^chr[0-9X]*\t' | cut -f1,2,3 > low_complexity.bed
grep 'Simple_repeat' repeatmask.txt  | grep $'^chr[0-9X]*\t' | cut -f1,2,3 > microsatellite.bed
grep 'DNA' repeatmask.txt  | grep $'^chr[0-9X]*\t' | cut -f1,2,3 > DNArep.bed #includes 'DNA?' ones, too

sort -k1,1 -k2,2n low_complexity.bed > low_complexity_sorted.bed 
sort -k1,1 -k2,2n DNArep.bed > DNArep_sorted.bed
sort -k1,1 -k2,2n microsatellite.bed > microsatellite_sorted.bed

rm low_complexity.bed
rm microsatellite.bed
rm DNArep.bed

bedtools merge -i microsatellite_sorted.bed -i low_complexity_sorted.bed -i DNArep_sorted.bed > region.bed

#running Snakefile
bsub < scripts/snakemake.sh
