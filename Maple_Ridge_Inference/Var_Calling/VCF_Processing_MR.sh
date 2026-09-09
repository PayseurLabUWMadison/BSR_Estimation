#!/bin/bash 


## Index the new vcf files for easier subsetting in the future steps 
../packages/gatk-4.2.0.0/gatk IndexFeatureFile -I maple_ridge_combined_sorted_merged_rmdups.vcf.gz



#Filter out variants that are close to genes using SelectVariants (several different buffers around the genes to explore how diversity changes as we get further from genic regions)
../packages/gatk-4.2.0.0/gatk SelectVariants -R ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -V maple_ridge_combined_sorted_merged_rmdups.vcf.gz -XL ../bed_files/merged_beds/P_man_merged_inv_repeats_10k_genes.bed --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs.vcf.gz


# Apply generic filters to each vcf with different levels of padding around genes 

../packages/gatk-4.2.0.0/gatk VariantFiltration -R ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna \
-V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs.vcf.gz  -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfilters.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 


## Remove the variants that do not have "PASS" in the provided filter sections 
../packages/gatk-4.2.0.0/gatk SelectVariants -R ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfilters.vcf.gz \
--exclude-filtered -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm.vcf.gz


##########
### Addtional Variant Processing 
##########


## Allele Depth Filtering (perl script from Peicheng Jing)
gunzip maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm.vcf.gz 
perl ../scripts/filter_for_AD.pl maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm.vcf maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass.vcf
# rezip and index 
bgzip -c maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass.vcf > maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass.vcf.gz
../packages/gatk-4.2.0.0/gatk IndexFeatureFile -I maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass.vcf.gz


## X Chromosome Maximum Depth Filtering 
# Process with seperate depth thresholds for males and females to account for differences in ploidy  
# Extract females from vcf
../packages/gatk-4.2.0.0/gatk SelectVariants -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass.vcf.gz -sn mr_females.list  -L NC_056031.1 -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Females.vcf.gz
# Extract males from vcf 
../packages/gatk-4.2.0.0/gatk SelectVariants -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass.vcf.gz -L NC_056031.1 -xl-sn mr_females.list -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Males.vcf.gz

## Depth Filter For Females (Using Average Depth To Set Threshold)
../packages/gatk-4.2.0.0/gatk VariantFiltration -R  ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Females.vcf.gz -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Females_DP.vcf.gz --genotype-filter-name "DP_Max46" --genotype-filter-expression "DP > 46"
 
## Depth Filter For Males (Using Average Depth To Set Threshold)
../packages/gatk-4.2.0.0/gatk VariantFiltration -R ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Males.vcf.gz -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Males_DP.vcf.gz --genotype-filter-name "DP_Max26" --genotype-filter-expression "DP > 26"

## Merge the filtered VCFs 
bcftools merge maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Females_DP.vcf.gz maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_Males_DP.vcf.gz -O z -o maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_All_DP.vcf.gz 
../packages/gatk-4.2.0.0/gatk IndexFeatureFile -I maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_All_DP.vcf.gz

## Remove variants that failed either male or female filter
../packages/gatk-4.2.0.0/gatk SelectVariants -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_All_DP.vcf.gz --max-filtered-genotypes 0 -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_All_maxDPrm.vcf.gz 

## Autosome Max Depth Filtering to be combined with 
../packages/gatk-4.2.0.0/gatk SelectVariants -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass.vcf.gz -L Autosomes.list -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes.vcf.gz
# Filter Autosomal variants for max depth (Using Average Depth To Set Threshold)
../packages/gatk-4.2.0.0/gatk VariantFiltration -R  ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes.vcf.gz -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_DP.vcf.gz --genotype-filter-name "DP_Max41" --genotype-filter-expression "DP > 41"
# Remove Variants that failed filter 
../packages/gatk-4.2.0.0/gatk SelectVariants -V maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_DP.vcf.gz --max-filtered-genotypes 0 -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm.vcf.gz 

## Combine X and Autosome VCFs back together 
../packages/gatk-4.2.0.0/gatk MergeVcfs -I maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_Autosomes_maxDPrm.vcf.gz -I maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_X_All_maxDPrm.vcf.gz -O maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_maxDPrm.vcf.gz



###################
### REMOVE RELATED INDIVIDUALS FROM VCFs 
################### 

# using relatednessfilter results from King and plinkQC, keep individuals that passed relatedness test for each population (list of samples comes from Filter_Related_GulfIslands.R run locally)

#### MAPLE RIDGE 
bcftools view maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_maxDPrm.vcf.gz -Oz -o unrelated_maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_maxDPrm.vcf.gz --samples "FSB14,FSB17,FSB2,FSB20,FSB26,FSB27,FSB3,FSB37,FSB4,FSB41,FSB42,FSB47,FSB5,FSB51,FSB52,FSB56,FSB7" 

../packages/gatk-4.2.0.0/gatk IndexFeatureFile -I unrelated_maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_maxDPrm.vcf.gz


## Seperate the vcfs for the X and autosomes for the next steps
#Maple Ridge
../packages/gatk-4.2.0.0/gatk SelectVariants -V unrelated_maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_maxDPrm.vcf.gz -R ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -L NC_056031.1 -O unrelated_maple_ridge_only_X.vcf.gz 
../packages/gatk-4.2.0.0/gatk SelectVariants -V unrelated_maple_ridge_combined_sorted_merged_rmdups_10kbgenes_bialSNPs_genericfiltersrm_ADPass_maxDPrm.vcf.gz -R ../reference_files/GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -XL NC_056031.1 -O unrelated_maple_ridge_only_Auto.vcf.gz 












