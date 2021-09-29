#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q long   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID
#BSUB -J snakemake
#BSUB -M 2048

module purge
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
module load python-3.9.5-gcc-9.3.0-jtayjft

conda init
. ~/.bashrc

conda activate Rabia
snakemake --latency-wait 300 --rerun-incomplete --cluster "bsub -M {resources.mem_mb} -n {threads} -q standard -J snakemake -e {params.err} -o {params.out}" --jobs 500 --keep-going --use-singularity
