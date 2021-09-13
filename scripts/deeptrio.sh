#!/usr/bin/env bash

#BSUB -n 49          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID
#BSUB -J deeptrio2X
#BSUB -M 81920

module purge
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp

cd ..

#try chr{{1..22},X}
singularity run docker://google/deepvariant:deeptrio-1.2.0  \
  /opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=WGS \
  --ref=GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --reads_child=crams/HG00408.final.cram \
  --reads_parent1=crams/HG00406.final.cram \
  --reads_parent2=crams/HG00407.final.cram \
  --output_vcf_child deeptrio_calls/trio2_chrX/HG00408.output.vcf.gz \
  --output_vcf_parent1 deeptrio_calls/trio2_chrX/HG00406.output.vcf.gz \
  --output_vcf_parent2 deeptrio_calls/trio2_chrX/HG00407.output.vcf.gz \
  --sample_name_child 'HG00408' \
  --sample_name_parent1 'HG00406' \
  --sample_name_parent2 'HG00407' \
  --num_shards 95  \
  --regions "chrX" \
  --intermediate_results_dir deeptrio_calls/trio2_chrX/intermediate_results_dir \
