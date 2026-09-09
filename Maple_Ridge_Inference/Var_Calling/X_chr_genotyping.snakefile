## Will Spurley
# 4/8/2025

# goal of this workflow: finish processing of variants following genotyping on CHTC (in parallel). Input gVCFs for each chromosome, output vcfs for the whole dataset 

import glob as glob
import os 


# Config File to Use:
configfile: "X_chr_genotyping.config.json"

# Input Functions:

def getGVCFs(pop, wildcard):
	inputs = []
	for sample in config[str(pop)]:
		inputs.append(config["GVCF_path"] + str(sample) + "_chr_" + str(wildcard) + "_sorted_merged_rmdups.g.vcf.gz")
	return inputs

def FormatGVCFs(pop, wildcard):
	inputs = []
	for sample in config[str(pop)]:
		inputs.append(config["GVCF_path"] + str(sample) + "_chr_" + str(wildcard) + "_sorted_merged_rmdups.g.vcf.gz")
	return " -V ".join(inputs)


def GatherVCFs(pop):
	inputs = []
	for chr in config["Chromosomes"]: 
		inputs.append(config["vcf_path"]+ str(pop)+ "_chr_" + str(chr) + "_sorted_merged_rmdups.vcf.gz")
	return " I=".join(inputs)


# -------------------------------------------------------------
# Rule all: determine what files to create using this workflow
# -------------------------------------------------------------

rule all: 
	input:
		expand((config["vcf_path"] + "maple_ridge_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config["Chromosomes"]),
		expand((config["vcf_path"] + "saturna_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config["Chromosomes"]),
		expand((config["vcf_path"] + "pender_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config["Chromosomes"]),
		expand((config["vcf_path"] + "maple_ridge_fem_only_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config["Chromosomes"]),
		(config["vcf_path"] + "maple_ridge_combined_sorted_merged_rmdups.vcf.gz"),
		(config["vcf_path"] + "saturna_combined_sorted_merged_rmdups.vcf.gz"),
		(config["vcf_path"] + "pender_combined_sorted_merged_rmdups.vcf.gz"),
		(config["vcf_path"] + "maple_ridge_fem_only_combined_sorted_merged_rmdups.vcf.gz")



# combine GVCFs for each population 
rule combine_GVCF_Maple_Ridge: 
	input:
		getGVCFs("Maple_Ridge", "{chrom}")
	output:
		outfile = (config["combined_gvcf_path"] + "maple_ridge_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["GVCF_path"],
		ref = config["ref_genome"],
		gvcfs = FormatGVCFs("Maple_Ridge", "{chrom}"),
		chrom = "{chrom}"
	shell: 
		"{params.gatk} CombineGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcfs} -L {params.chrom} -O {output.outfile}"


rule combine_GVCF_Saturna: 
	input:
		getGVCFs("Saturna", "{chrom}")
	output:
		outfile = (config["combined_gvcf_path"] + "saturna_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["GVCF_path"],
		ref = config["ref_genome"],
		gvcfs = FormatGVCFs("Saturna", "{chrom}"),
		chrom = "{chrom}"
	shell: 
		"{params.gatk} CombineGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcfs} -L {params.chrom} -O {output.outfile}"


rule combine_GVCF_Pender: 
	input:
		getGVCFs("Pender", "{chrom}")
	output:
		outfile = (config["combined_gvcf_path"] + "pender_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["GVCF_path"],
		ref = config["ref_genome"],
		gvcfs = FormatGVCFs("Pender", "{chrom}"),
		chrom = "{chrom}"
	shell: 
		"{params.gatk} CombineGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcfs} -L {params.chrom} -O {output.outfile}"


rule combine_GVCF_MR_Females:
	input:
		getGVCFs("MR_Females", "{chrom}")
	output:
		outfile = (config["combined_gvcf_path"] + "maple_ridge_fem_only_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["GVCF_path"],
		ref = config["ref_genome"],
		gvcfs = FormatGVCFs("MR_Females", "{chrom}"),
		chrom = "{chrom}"
	shell: 
		"{params.gatk} CombineGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcfs} -L {params.chrom} -O {output.outfile}"



## Joint variant calling using GVCFs
rule genotype_Maple_Ridge:
	input:
		(config["combined_gvcf_path"] + "maple_ridge_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	output:
		(config["vcf_path"] + "maple_ridge_chr_{chrom}_sorted_merged_rmdups.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["combined_gvcf_path"],
		ref = config["ref_genome"],
		gvcf = (config["combined_gvcf_path"] + "maple_ridge_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz"),
		outfile = (config["vcf_path"] + "maple_ridge_chr_{chrom}_sorted_merged_rmdups.vcf.gz"),
		chrom = "{chrom}"
	shell:
		"{params.gatk} GenotypeGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcf} -L {params.chrom} -O {params.outfile}"


rule genotype_Pender:
	input:
		(config["combined_gvcf_path"] + "pender_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	output:
		(config["vcf_path"] + "pender_chr_{chrom}_sorted_merged_rmdups.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["combined_gvcf_path"],
		ref = config["ref_genome"],
		gvcf = (config["combined_gvcf_path"] + "pender_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz"),
		outfile = (config["vcf_path"] + "pender_chr_{chrom}_sorted_merged_rmdups.vcf.gz"),
		chrom = "{chrom}"
	shell:
		"{params.gatk} GenotypeGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcf} -L {params.chrom} -O {params.outfile}"


rule genotype_Saturna:
	input:
		(config["combined_gvcf_path"] + "saturna_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	output:
		(config["vcf_path"] + "saturna_chr_{chrom}_sorted_merged_rmdups.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["combined_gvcf_path"],
		ref = config["ref_genome"],
		gvcf = (config["combined_gvcf_path"] + "saturna_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz"),
		outfile = (config["vcf_path"] + "saturna_chr_{chrom}_sorted_merged_rmdups.vcf.gz"),
		chrom = "{chrom}"
	shell:
		"{params.gatk} GenotypeGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcf} -L {params.chrom} -O {params.outfile}"


rule genotype_MR_Females:
	input:
		(config["combined_gvcf_path"] + "maple_ridge_fem_only_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz")
	output:
		(config["vcf_path"] + "maple_ridge_fem_only_chr_{chrom}_sorted_merged_rmdups.vcf.gz")
	params:
		gatk = config["gatk_path"],
		tmp_dir = config["combined_gvcf_path"],
		ref = config["ref_genome"],
		gvcf = (config["combined_gvcf_path"] + "maple_ridge_fem_only_chr_{chrom}_sorted_merged_rmdups.g.vcf.gz"),
		outfile = (config["vcf_path"] + "maple_ridge_fem_only_chr_{chrom}_sorted_merged_rmdups.vcf.gz"),
		chrom = "{chrom}"
	shell:
		"{params.gatk} GenotypeGVCFs --tmp-dir {params.tmp_dir} -R {params.ref} -V {params.gvcf} -L {params.chrom} -O {params.outfile}"


## Combine all chromosomal VCFs into population specific VCFs 
rule gather_MR:
	input:
		expand((config["vcf_path"] + "maple_ridge_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config['Chromosomes'])
	output:
		(config["vcf_path"] + "maple_ridge_combined_sorted_merged_rmdups.vcf.gz")
	params:
		vcfs = GatherVCFs("maple_ridge"),
		outfile = (config["vcf_path"]+"maple_ridge_combined_sorted_merged_rmdups.vcf.gz"),
		gatk = config['gatk_path']
	shell:
		"{params.gatk} GatherVcfs I={params.vcfs} O={params.outfile}"


rule gather_Pender:
	input:
		expand((config["vcf_path"] + "pender_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config["Chromosomes"])
	output:
		(config["vcf_path"] + "pender_combined_sorted_merged_rmdups.vcf.gz")
	params:
		vcfs = GatherVCFs("pender"),
		outfile = (config["vcf_path"]+"pender_combined_sorted_merged_rmdups.vcf.gz"),
		gatk = config['gatk_path']
	shell:
		"{params.gatk} GatherVcfs I={params.vcfs} O={params.outfile}"

rule gather_Saturna:
	input:
		expand((config["vcf_path"] + "saturna_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config['Chromosomes'])
	output:
		(config["vcf_path"] + "saturna_combined_sorted_merged_rmdups.vcf.gz")
	params:
		vcfs = GatherVCFs("saturna"),
		outfile = (config["vcf_path"]+"saturna_combined_sorted_merged_rmdups.vcf.gz"),
		gatk = config['gatk_path']
	shell:
		"{params.gatk} GatherVcfs I={params.vcfs} O={params.outfile}"

rule gather_MR_Females:
	input:
		expand((config["vcf_path"] + "maple_ridge_fem_only_chr_{chrom}_sorted_merged_rmdups.vcf.gz"), chrom=config["Chromosomes"])
	output:
		(config["vcf_path"] + "maple_ridge_fem_only_combined_sorted_merged_rmdups.vcf.gz")
	params:
		vcfs = GatherVCFs("maple_ridge_fem_only"),
		outfile = (config["vcf_path"]+"maple_ridge_fem_only_combined_sorted_merged_rmdups.vcf.gz"),
		gatk = config['gatk_path']
	shell:
		"{params.gatk} GatherVcfs I={params.vcfs} O={params.outfile}"









