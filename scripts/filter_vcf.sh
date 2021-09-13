#!/usr/bin/env bash

#BSUB -n 1          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e /dev/null   # error file name in which %J is replaced by the job ID
#BSUB -o /dev/null  # output file name in which %J is replaced by the job ID
#BSUB -J vcf_filter_ms


cd /hps/nobackup/goldman/denovo_switch/octopus_calls

#while read i; do bcftools view -H -e 'INFO/DENOVO!=1 | FILTER!="PASS" | FMT/FT!="PASS" | FMT/DP <30' -T ^../low_complexity.bed trio${i}.vcf.gz | wc -l >> denovo/filtered/DP30_wo_lc.txt; done < ../logs/octopus/successful.txt

#while read i; do bcftools view -e 'INFO/DENOVO!=1 | FILTER!="PASS" | FMT/FT!="PASS" | FMT/DP <24' -T ^../low_complexity.bed trio${i}.vcf.gz | bcftools view -H -T ^../microsatellite.bed | wc -l >>denovo/filtered/DP24_wo_lc_ms.txt; done < ../logs/octopus/successful.txt

#while read i; do bcftools view -e 'INFO/DENOVO!=1 | FILTER!="PASS" | FMT/FT!="PASS" | FMT/DP <24' -T ^../low_complexity.bed trio${i}.vcf.gz | bcftools view -T ^../microsatellite.bed | bcftools view -H -T ^../DNArep.bed |wc -l >>denovo/filtered/DP24_wo_lc_ms_dr.txt; done < ../logs/octopus/successful.txt

#for i in 225 147 277 269; do bcftools view -e 'INFO/DENOVO!=1 | FILTER!="PASS" | FMT/FT!="PASS" | FMT/DP <24' -T ^../low_complexity.bed trio${i}.vcf.gz | bcftools view -T ^../microsatellite.bed | bcftools view -T ^../DNArep.bed >denovo/filtered/trio${i}_DP24_wo_lc_ms_dr.vcf; done
#for i in 81 287; do bcftools view -e 'INFO/DENOVO!=1 | FILTER!="PASS" | FMT/FT!="PASS" | FMT/DP <24' -T ^../low_complexity.bed trio${i}.vcf.gz | bcftools view -T ^../microsatellite.bed | bcftools view -T ^../DNArep.bed >denovo/filtered/trio${i}_DP24_wo_lc_ms_dr.vcf; done

while read i; do bcftools view -H -e 'INFO/DENOVO!=1 | FILTER!="PASS" | FMT/FT!="PASS" | FMT/DP <20 | MP<40 ' -T ^../region.bed trio${i}.vcf.gz |wc -l >>denovo/filtered/DP20_MP40_wo_lc_ms_dr.txt; done < ../logs/octopus/successful.txt

