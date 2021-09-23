configfile: "config.yaml"

localrules: all

wildcard_constraints:
    trio="\d+"

#seq=range(1,603)  #trio index
#seq=range(1,301) #Rabia
#seq=range(301,603) #Conor


#trios=[x for x in seq]
trios=269

rule all:
	input:
		#expand("octopus_calls/trio{trio}.vcf.gz",trio=trios),
		expand("deeptrio_calls/trio{trio}_denovo.vcf.gz",trio=trios),
		expand("denovo_variants/octopus/trio{trio}.vcf",trio=trios)

rule samtools_index:
	input:
		#sample=lambda wildcards: "crams/"+wildcards.sample+".cram"
		sample="crams/{sample}.cram"


	output:
		"crams/{sample}.cram.crai"

	params:
		err=lambda wildcards: "logs/samtools_index/err." +wildcards.sample,
		out=lambda wildcards: "logs/samtools_index/out." +wildcards.sample

	threads: 2

	resources:
		mem_mb=2048

	shell:
		"samtools index {input.sample}"


rule octopus:
	input:
		"GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
		midx=lambda wildcards: "crams/"+ config['trios'][wildcards.trio]["mother"]+".final.cram.crai",
		fidx=lambda wildcards: "crams/"+ config['trios'][wildcards.trio]["father"]+".final.cram.crai",
		cidx=lambda wildcards: "crams/"+ config['trios'][wildcards.trio]["child"]+".final.cram.crai",
		ref="GRCh38_full_analysis_set_plus_decoy_hla.fa",
		mother=lambda wildcards: "crams/" +config['trios'][wildcards.trio]["mother"]+".final.cram",
		father=lambda wildcards: "crams/" +config['trios'][wildcards.trio]["father"]+".final.cram",
		child=lambda wildcards: "crams/" +config['trios'][wildcards.trio]["child"]+".final.cram"


	output:
		"octopus_calls/trio{trio}.vcf.gz"

	params:
		mum=lambda wildcards: config['trios'][wildcards.trio]["mother"],
		dad=lambda wildcards: config['trios'][wildcards.trio]["father"],
		err=lambda wildcards: "logs/octopus/err.trio" +wildcards.trio,
		out=lambda wildcards: "logs/octopus/out.trio" +wildcards.trio
	
	threads: 95

	resources:
		mem_mb=91200

	shell:
		"octopus -R {input.ref} -I {input.father} {input.mother} {input.child} -M {params.mum} -F {params.dad} -o {output} --threads {threads} \		--sequence-error-model PCR-FREE.NOVASEQ --read-linkage PAIRED -T chr1 to chrX"

rule octopusFilter:
	input:
		vcf="octopus_calls/trio{trio}.vcf.gz",
		region="region.bed"

	output:
		"denovo_variants/octopus/trio{trio}.vcf"

	params:
		err=lambda wildcards: "logs/filter/err.trio" +wildcards.trio,
		out=lambda wildcards: "logs/filter/out.trio" +wildcards.trio
		#err="logs",
		#out="/dev/null"

	threads: 1

	resources:
		mem_mb=2000

	shell:
		"bcftools view -e 'INFO/DENOVO!=1 | INFO/REVERSION==1 | FILTER!=\"PASS\" | FMT/FT!=\"PASS\" | FMT/DP <24 | MP<40 ' -T ^{input.region} {input.vcf} > {output}"

	
rule deepTrio:
	input:
		"GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
		midx=lambda wildcards: "crams/"+ config['trios'][wildcards.trio]["mother"]+".final.cram.crai",
		fidx=lambda wildcards: "crams/"+ config['trios'][wildcards.trio]["father"]+".final.cram.crai",
		cidx=lambda wildcards: "crams/"+ config['trios'][wildcards.trio]["child"]+".final.cram.crai",
		ref="GRCh38_full_analysis_set_plus_decoy_hla.fa",
		mother=lambda wildcards: "crams/" +config['trios'][wildcards.trio]["mother"]+".final.cram",
		father=lambda wildcards: "crams/" +config['trios'][wildcards.trio]["father"]+".final.cram",
		child=lambda wildcards: "crams/" +config['trios'][wildcards.trio]["child"]+".final.cram"
	output:
		mothervcf="deeptrio_calls/trio{trio}/mum.output.vcf.gz",
		fathervcf="deeptrio_calls/trio{trio}/dad.output.vcf.gz",
		childvcf="deeptrio_calls/trio{trio}/child.output.vcf.gz",
		mothergvcf="deeptrio_calls/trio{trio}/mum.g.vcf.gz",
		fathergvcf="deeptrio_calls/trio{trio}/dad.g.vcf.gz",
		childgvcf="deeptrio_calls/trio{trio}/child.g.vcf.gz",
		int=directory("deeptrio_calls/trio{trio}/intermediate_results_dir")
	params:
		mum=lambda wildcards: config['trios'][wildcards.trio]["mother"],
		dad=lambda wildcards: config['trios'][wildcards.trio]["father"],
		child=lambda wildcards: config['trios'][wildcards.trio]["child"],
		err=lambda wildcards: "logs/deeptrio/err.trio" +wildcards.trio,
		out=lambda wildcards: "logs/deeptrio/out.trio" +wildcards.trio

	threads: 90

	resources:
		mem_mb=81920

	container:
		"docker://google/deepvariant:deeptrio-1.2.0"

	shell:
		"/opt/deepvariant/bin/deeptrio/run_deeptrio \
		--model_type=WGS \
		--ref={input.ref} \
		--reads_child={input.child} \
		--reads_parent1={input.mother} \
		--reads_parent2={input.father} \
		--output_vcf_child {output.childvcf} \
		--output_vcf_parent1 {output.mothervcf} \
		--output_vcf_parent2 {output.fathervcf} \
		--sample_name_child '{params.child}' \
		--sample_name_parent1 '{params.mum}' \
		--sample_name_parent2 '{params.dad}' \
		--num_shards {threads}  \
		--regions \"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
		chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22\" \
		--intermediate_results_dir {output.int} \
		--output_gvcf_child {output.childgvcf} \
		--output_gvcf_parent1 {output.mothergvcf} \
		--output_gvcf_parent2 {output.fathergvcf}"


rule glnexus_merge:
	input:
		m=lambda wildcards: "deeptrio_calls/trio"+ wildcards.trio+"/mum.g.vcf.gz",
		f=lambda wildcards: "deeptrio_calls/trio"+ wildcards.trio+"/dad.g.vcf.gz",
		c=lambda wildcards: "deeptrio_calls/trio"+ wildcards.trio+"/child.g.vcf.gz"
	output:
		"deeptrio_calls/trio{trio}.vcf.gz"
	container:
		"docker://quay.io/mlin/glnexus:v1.3.1"
	threads: 10
	resources:
		mem_mb=10000
	params:
		err=lambda wildcards: "logs/glnexus/err.trio" +wildcards.trio,
		out=lambda wildcards: "logs/glnexus/out.trio" +wildcards.trio

	shell:
		"/usr/local/bin/glnexus_cli -t {threads} -m {resources.mem_mb} \
		--config DeepVariant_unfiltered {input.f} {input.m} {input.c} | bcftools view |bgzip -c > {output}"


rule DeepTrioDenovo:
	input:
		"deeptrio_calls/trio{trio}.vcf.gz"
	output:
		"deeptrio_calls/trio{trio}_denovo.vcf.gz"

	threads: 1
		
	resources:
		mem_mb=2000
	params:
		child=lambda wildcards: config['trios'][wildcards.trio]["child"],
		err=lambda wildcards: "logs/DTdenovo/err.trio" +wildcards.trio,
		out=lambda wildcards: "logs/DTdenovo/out.trio" +wildcards.trio
	shell:
		"python scripts/denovo.py all all {input} {params.child} {output}"
