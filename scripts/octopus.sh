#!/usr/bin/env bash

#BSUB -n 95          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID
#BSUB -J octopus_autosomes_X_slow
#BSUB -M 71680

cd ..
conda init
. ~/.bashrc
conda activate Rabia
octopus -R GRCh38_full_analysis_set_plus_decoy_hla.fa -I crams/HG00405.final.cram crams/HG00403.final.cram crams/HG00404.final.cram -M HG00404 -F HG00403 -o octopus_calls/trio1_autosomes_X_slow.vcf.gz --sequence-error-model PCR-FREE.NOVASEQ --read-linkage PAIRED --threads 95 -T chr1 to chrX
