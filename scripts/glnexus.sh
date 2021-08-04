#!/usr/bin/env bash

#BSUB -n 95          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID
#BSUB -J deeptrio_autosomes_X
#BSUB -M 81920

module purge
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp

cd ..

singularity run deepvariant_deeptrio-1.1.0.sif  \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=WGS \
  --ref=GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --reads_child=crams/HG00405.final.cram \
  --reads_parent1=crams/HG00403.final.cram \
  --reads_parent2=crams/HG00404.final.cram \
  --output_vcf_child deeptrio_calls/trio1/HG00405.output.vcf.gz \
  --output_vcf_parent1 deeptrio_calls/trio1/HG00403.output.vcf.gz \
  --output_vcf_parent2 deeptrio_calls/trio1/HG00404.output.vcf.gz \
  --sample_name_child 'HG00405' \
  --sample_name_parent1 'HG00403' \
  --sample_name_parent2 'HG00404' \
  --num_shards $(nproc)  \
  --regions "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX" \
  --intermediate_results_dir deeptrio_calls/trio1/intermediate_results_dir \
  --output_gvcf_child deeptrio_calls/trio1/HG00405.g.vcf.gz \
  --output_gvcf_parent1 deeptrio_calls/trio1/HG00403.g.vcf.gz \
  --output_gvcf_parent2 deeptrio_calls/trio1/HG00404.g.vcf.gz
