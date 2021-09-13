#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID
#BSUB -J pipeline

#Data: High coverage 1KG trios
#https://www.internationalgenome.org/data-portal/data-collection/30x-grch38

#create a directory for data manipulation
if [ ! -d data ]; then
    mkdir data;
fi

cd data

#fetch pedigree file
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

#data manipulation
ml load r-4.0.3-gcc-9.3.0-4l6eluj

rwd=$(pwd) # r working directory
Rscript ../scripts/data.R $rwd

#create crams directory to store cram files
if [ ! -d crams ]; then
    mkdir crams;
fi

#download the table containing the cram file links from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
#You can also get the same list from the github repo:
wget https://github.com/rabiafidan/denovoTS/blob/b029b11f4badec37d1493ac570b98ae567fb05c9/igsr_30x%20GRCh38.tsv

#obtain cram file links of the individuals of the trios
grep cram igsr_30x\ GRCh38.tsv | cut -f1 >links.txt
for i in $(cat inds.txt); do grep $i links.txt >> links2.txt; done
mv links2.txt links.txt

cd ../scripts  #MANAGE PATHS AND GIVE IT AS ARGUMENTS TO DOWNLOAD CODE
#download cram files and the reference fasta
bsub -J "dwnld[1-1793]%40" < download.sh

#index the reference fasta
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

#create the configuration file for Snakefile
ml load python-3.9.0-gcc-9.3.0-5t75egs
python 


