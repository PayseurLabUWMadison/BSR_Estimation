## Will Spurley
# created on: 2/28/2025

# Goal of this script: Create a pipeline for processing raw short read data into bam files ready for genotyping, specifically with the intent of comparing the X chromosome to the rest of the genome

# import libraries (snakemake works in python language)
import glob as glob
import os

# Config File to Use:
configfile: 'X_chr_ngs_workflow_batch3.config.json'


#------------------------------
# Input Functions and rule all
#------------------------------

# The goal of functions is to quickly build a list of input files that can be parsed through their respective rules to generate the correct output files
	# example: grab all fastq files in the fastq directory to generate 

# create various functions as needed 
def getFastQCInputsR1(wildcard): # Rework to use two wildcards, making the S### part of the file a * 
	return [os.path.basename(x) for x in glob.glob((str(config["fastq_path"]) + str(wildcard) + "_R1_001" + config["fastq_suffix"]))]
def getFastQCInputsR2(wildcard):
	return [os.path.basename(x) for x in glob.glob((str(config["fastq_path"]) + str(wildcard) + "_R2_001" + config["fastq_suffix"]))]

def getTrimInputsR1(wildcard):
	return [os.path.basename(x) for x in glob.glob((str(config["fastq_path"]) + str(wildcard[0]) + "_*_" + str(wildcard[1]) + "_R1_001" + config["fastq_suffix"]))]
def getTrimInputsR2(wildcard):
	return [os.path.basename(x) for x in glob.glob((str(config["fastq_path"]) + str(wildcard[0]) + "_*_" + str(wildcard[1]) + "_R2_001" + config["fastq_suffix"]))]


def getFastQCTargets():
	r1_targets = list()
	r2_targets = list()
	for sample in config["sample"]:
		for lane in config["lane_id_list"]:
			r1_targets.extend([(str(config["fastqc_report_path"]) + os.path.basename(x).split(".")[0] + "_fastqc.zip") for x in glob.glob((str(config["fastq_path"]) + sample + "_*_" + lane + "_R1_001" + str(config["fastq_suffix"])))])
			r2_targets.extend([(str(config["fastqc_report_path"]) + os.path.basename(x).split(".")[0] + "_fastqc.zip") for x in glob.glob((str(config["fastq_path"]) + sample + "_*_" + lane + "_R2_001" + str(config["fastq_suffix"])))])
	return r1_targets, r2_targets


def getSamtoolsInput(wildcard):
	inputs = []
	for lane in config["lane_id_list"]:
		inputs.append(config['bam_path'] + str(wildcard) + "_" + lane + "_sorted.bam")
	return inputs


# rule all: designate the desired outputs (what are the end goals and what rules does snakemake need to get to that point?)
rule all:
	input:
		#initial Fast QC reports 
		getFastQCTargets(),
		#post-adapter trim fastQC reports
		expand((config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_1P_fastqc.zip"), sample=config["sample"], lane_id_list=config["lane_id_list"]),
		expand((config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_2P_fastqc.zip"), sample=config["sample"], lane_id_list=config["lane_id_list"]),
		# # post-quality trim fastQC reports 
		expand((config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_1P_fastqc.zip"), sample=config["sample"], lane_id_list=config["lane_id_list"]),
		expand((config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_2P_fastqc.zip"), sample=config["sample"], lane_id_list=config["lane_id_list"]),
		# Summary statistics on intial alignment 
		expand((config["samtools_report_path"] + "{sample}_{lane_id_list}_sorted.flagstat.txt"), sample=config["sample"], lane_id_list=config["lane_id_list"]),
		# Mark duplicates report 
		expand((config["picard_report_path"] + "{sample}_sorted_merged_DuplicationMetrics.txt"), sample=config["sample"]),
		# XYalign output for read depth
		#expand((config["xyalign_plot_path"] + "{sample}_noprocessing_{chromosome}_Depth_GenomicScatter.pdf"), sample=config["sample"], chromosome=config["xyalign_chrs"])


#----------------------------------
# QC steps to prepare for alignment
#----------------------------------


# intial FastQC for comparison with downstream steps 
rule vis_fastqc_init:
	input:
	output:
		(config["fastqc_report_path"] + "{wildcard}_R1_001_" + "fastqc.zip"),
		(config["fastqc_report_path"] + "{wildcard}_R2_001_" + "fastqc.zip")
	params:
		indir = config["fastq_path"],
		wkdir = config["working_path"],
		fastq_r1 = getFastQCInputsR1, 
		fastq_r2 = getFastQCInputsR2,
		outdir = config["fastqc_report_path"],
		fastqc = config["fastqc_script"]
	shell:
		"""
		cd "{params.indir}"
		{params.fastqc} {params.fastq_r1} {params.fastq_r2} --outdir={params.outdir}
		cd "{params.wkdir}"
		"""

## Trim the adapter sequences off the reads using bbduk 
rule trim_adapters: 
	input: 
	output:
		(config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_1P.fastq.gz"),
		(config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_2P.fastq.gz")
	params:
		fastq_r1 = getTrimInputsR1,
		fastq_r2 = getTrimInputsR2,
		bbduk = config["bbduk_script"],
		indir = config["fastq_path"],
		wkdir = config["working_path"],
		adapter = config["adapter"],
		fastq_r1_out = str(config["filtered_fastq_path"]) + "{sample}_{lane_id_list}_adaptrim_1P.fastq.gz",
		fastq_r2_out = str(config["filtered_fastq_path"]) + "{sample}_{lane_id_list}_adaptrim_2P.fastq.gz"
	shell:
		"""
		cd "{params.indir}"
		{params.bbduk} -Xmx1g in1={params.fastq_r1} in2={params.fastq_r2} out1={params.fastq_r1_out} out2={params.fastq_r2_out} ref={params.adapter} ktrim=r k=23 mink=11 hdist=1 tpe tbo
		cd {params.wkdir}		
		"""


## Perform another round of QC to examine how trimming adapters changes quality 
rule vis_fastqc_trim:
	input:
		fastq_r1 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_1P.fastq.gz"),
		fastq_r2 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_2P.fastq.gz")
	output:
		(config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_1P_fastqc.zip"),
		(config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_2P_fastqc.zip")
	params:
		outdir = config["fastqc_report_path"],
		fastqc = config["fastqc_script"]
	shell:
		"{params.fastqc} {input.fastq_r1} {input.fastq_r2} --outdir={params.outdir}"


## Trimmomatic to remove low quality reads on ends
rule qual_trim: 
	input:
		fastq_r1 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_1P.fastq.gz"),
		fastq_r2 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_2P.fastq.gz")
	output:
		(config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_1P.fastq.gz"),
		(config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_1U.fastq.gz"),
		(config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_2P.fastq.gz"),
		(config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_2U.fastq.gz")
	params:
		trimmomatic = config["trimmomatic_script"],
		baseout = config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim.fastq.gz"
	shell:
		"java -jar {params.trimmomatic} PE {input.fastq_r1} {input.fastq_r2} -baseout {params.baseout} SLIDINGWINDOW:4:15 MINLEN:36"


## FastQC on the adapter- and quality-trimmed reads (goal for this is to see an increase in first Fastqc stat)
rule vis_fastqc_qualtrim:
	input:
		fastq_r1 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_1P.fastq.gz"),
		fastq_r2 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_2P.fastq.gz")
	output:
		(config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_1P_fastqc.zip"),
		(config["fastqc_report_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_2P_fastqc.zip")
	params:
		fastqc = config["fastqc_script"],
		outdir = config["fastqc_report_path"]
	shell:
		"{params.fastqc} {input.fastq_r1} {input.fastq_r2} --outdir={params.outdir}"


#------------------
# Alignment and QC
#------------------

# Intial alignment to the reference for each lane 
rule init_alignment: 
	input:
		fastq_r1 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_1P.fastq.gz"),
		fastq_r2 = (config["filtered_fastq_path"] + "{sample}_{lane_id_list}_adaptrim_qualtrim_2P.fastq.gz")
	output:
		(config["bam_path"] + "{sample}_{lane_id_list}_sorted.bam")
	params:
		read_groups = config["read_groups_script"],
		ref_genome = config["ref_genome"],
		sample = "{sample}",
		out_path = (config["bam_path"] + "{sample}_{lane_id_list}_sorted.bam")
	shell:
		"bash {params.read_groups} {params.ref_genome} {input.fastq_r1} {input.fastq_r2} {params.sample} {params.out_path}"

rule samtool_stats:
	input:
		bam=(config["bam_path"] + "{sample}_{lane_id_list}_sorted.bam")
	output:
		(config["samtools_report_path"] + "{sample}_{lane_id_list}_sorted.flagstat.txt")
	params:
		output = (config["samtools_report_path"] + "{sample}_{lane_id_list}_sorted.flagstat.txt")
	shell:
		"samtools flagstat {input.bam} > {params.output}" 

rule sample_merge: 
	input:
		getSamtoolsInput
	output:
		(config["bam_path"] + "{sample}_sorted_merged.bam")
	params:
		output = (config["bam_path"] + "{sample}_sorted_merged.bam")
	shell:
		"samtools merge {params.output} {input}"


rule mark_duplicates:
	input:
		merged_bam = (config["bam_path"] + "{sample}_sorted_merged.bam")
	output:
		(config["bam_path"] + "{sample}_sorted_merged_rmdups.bam"),
		(config["picard_report_path"] + "{sample}_sorted_merged_DuplicationMetrics.txt")
	params:
		gatk = config["gatk_path"],
		new_bam = (config["bam_path"] + "{sample}_sorted_merged_rmdups.bam"),
		report = (config["picard_report_path"] + "{sample}_sorted_merged_DuplicationMetrics.txt"),
		tmp_dir = config["bam_path"],
		to_rm = (config["bam_path"] + "{sample}_sorted_merged.bam")
	shell:
		"""
		set -euxo pipefail
		{params.gatk} MarkDuplicates --TMP_DIR {params.tmp_dir} -I {input.merged_bam} -O {params.new_bam} -M {params.report} --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 170
		if [ -s {params.new_bam} ]
		then 
			rm {params.to_rm}
		fi 
		"""

rule compare_X_chrom:
	input:
		in_bam = (config["bam_path"] + "{sample}_sorted_merged_rmdups.bam")
	output:
		(config["xyalign_plot_path"] + "{sample}_noprocessing_{chromosome}_Depth_GenomicScatter.pdf")
	conda:
		"xyalign_env"
	params:
		ref_genome = config["ref_genome"],
		xyalign_out_path = config["xyalign_output_path"],
		X_chr = config["X_chr"],
		Autosome = config["Autosome_for_comp"]
	shell:
		"xyalign --CHARACTERIZE_SEX_CHROMS --ref {params.ref_genome} --bam {input.in_bam} --output_dir {params.xyalign_out_path} --sample_id {wildcards.sample} --cpus 1 --window_size 10000 --chromosomes {params.X_chr} {params.Autosome}  --x_chromosome {params.X_chr}"




## HaplotypeCaller must be performed in CHTC to speed the process up. Make sure this script is doing everything required before that 


