#!/bin/bash

set -e

# # Rename Chromosomes to make the downstream process smoother
bcftools annotate --rename-chrs P_man_rename.map -Oz -o maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm_rename.vcf.gz maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm.vcf.gz


## Convert Filtered VCFs to binary files for use in plink
#	Maple Ridge
plink1.9 --vcf maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm_rename.vcf.gz --make-bed --chr-set 22 no-xy --out maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm_dup

## Create LD pruned Autosomal variant callset for each population to identify related individuals 

# 	Maple Ridge
plink1.9 --bfile maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm_dup --indep-pairwise 50 10 0.1 --out maple_ridge_Auto_LDPrune_dup
plink1.9 --bfile maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm_dup --extract maple_ridge_Auto_LDPrune.prune.in --make-bed --out maple_ridge_Auto_LDPrune_dup

## run pairwise kinship analysis for each population 
#	Maple Ridge  
../packages/king -b maple_ridge_Auto_LDPrune_dup.bed --kinship --prefix maple_ridge_dup
