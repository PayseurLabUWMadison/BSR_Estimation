#!/bin/bash

set -e

IND=$1
CHR=$2
SEX=$3

BAM=${IND}_${CHR}_sorted_merged_rmdups.bam
BAI=${IND}_${CHR}_sorted_merged_rmdups.bam.bai
REF=GCF_003704035.1_HU_Pman_2.1.3_genomic.fna


#extract gatk from tarball
tar -xzf gatk-4.2.0.0.tar.gz

#Allow for java to run
tar -xzf openjdk-24_linux-x64_bin.tar.gz
export JAVA_HOME=$PWD/jdk-24
export PATH=$JAVA_HOME/bin:$PATH


if [[ "${CHR}" == "NC_056031.1" ]] && [[ "${SEX}" == "M" ]]; then
	#statements
	echo "Calling Variants as Haploid using: "
	echo " gatk-4.2.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller -R ${REF} -I ${BAM} -L ${CHR} -ploidy 1 -ERC GVCF -O ${IND}_chr_${CHR}_sorted_merged_rmdups.g.vcf.gz"
	# Variant calling for haploid chromosomes (Only X for males)
	gatk-4.2.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller -R ${REF} -I ${BAM} -L ${CHR} -ploidy 1 -ERC GVCF -O ${IND}_chr_${CHR}_sorted_merged_rmdups.g.vcf.gz

else
	echo "Calling Variants as Diploid using: "
	echo "gatk-4.2.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller -R ${REF} -I ${BAM} -L ${CHR} -ploidy 2 -ERC GVCF -O ${IND}_chr_${CHR}_sorted_merged_rmdups.g.vcf.gz"
	# Variant calling for diploid chromosomes (i.e. all autosomes and female X chromosomes)
	gatk-4.2.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller -R ${REF} -I ${BAM} -L ${CHR} -ploidy 2 -ERC GVCF -O ${IND}_chr_${CHR}_sorted_merged_rmdups.g.vcf.gz
fi


