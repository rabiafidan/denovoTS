#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q long   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID
#BSUB -J s_failed
#BSUB -M 2048


cd ..

conda init
. ~/.bashrc

conda activate Rabia
snakemake --cluster "bsub -M {resources.mem_mb} -n {threads} -q standard -J snakemake -e {params.err} -o {params.out}" --jobs 200
