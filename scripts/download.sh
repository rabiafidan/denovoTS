#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e /dev/null   # error file name in which %J is replaced by the job ID
#BSUB -o /dev/null  # output file name in which %J is replaced by the job ID



###retreive the individual link and download
llink=$(awk -v idx="${LSB_JOBINDEX}" 'FNR==idx {print $0}' links.txt)

cd ../crams/
wget $llink


