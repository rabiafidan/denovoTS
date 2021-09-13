configfile: "config.yaml"

#seq=range(1,603)  #trio index
#seq=range(2,301) #Rabia
#seq=range(301,603) #Conor


#trios=[x for x in seq]
trios=[364,370,390,413,432,511,513,530,546,595]

rule all:
	input:
		expand("octopus_calls/trio{trio}.vcf.gz",trio=trios)



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

	#wrapper:
        	#"0.77.0/bio/samtools/index"


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
                mum=lambda wildcards, input: input["mother"][6:-11],
                dad=lambda wildcards, input: input["father"][6:-11],
		err=lambda wildcards: "logs/octopus/err.trio" +wildcards.trio,
		out=lambda wildcards: "logs/octopus/out.trio" +wildcards.trio
	
	threads: 95

	resources:
		mem_mb=91200

	shell:
                "octopus -R {input.ref} -I {input.father} {input.mother} {input.child} -M {params.mum} -F {params.dad} -o {output} --threads {threads} \		--sequence-error-model PCR-FREE.NOVASEQ --read-linkage PAIRED -T chr1 to chrX"
