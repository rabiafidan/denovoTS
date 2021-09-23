#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e /dev/null   # error file name in which %J is replaced by the job ID
#BSUB -o /dev/null  # output file name in which %J is replaced by the job ID
#BSUB -J vcf_filter_ms


cd /hps/nobackup/goldman/denovo_switch/octopus_calls

for i in {1..602}; do ch=$(awk -v line="$i" 'FNR==line {print$0}' ../data/children.txt); bcftools view -e 'INFO/DENOVO!=1 | INFO/REVERSION==1 | FILTER!="PASS" | FMT/FT!="PASS" | FMT/DP <24 | MP<40 ' -T ^../region.bed trio${i}.vcf.gz >../denovo_variants/trio${i}_${ch}.vcf; done

