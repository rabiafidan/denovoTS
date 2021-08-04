#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e missing_link_error.%J.%I   # error file name in which %J is replaced by the job ID
#BSUB -o missing_link_output.%J.%I  # output file name in which %J is replaced by the job ID



###inds
llink=$(awk -v idx="${LSB_JOBINDEX}" 'FNR==idx {print $0}' /hps/nobackup/goldman/denovo_switch/links.txt)

cd /hps/nobackup/goldman/denovo_switch/crams/
wget $llink


###ref
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
