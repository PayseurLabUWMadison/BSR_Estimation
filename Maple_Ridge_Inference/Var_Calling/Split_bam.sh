#!/bin/bash


set -e


IND=$1
CHR=$2
BAM=${IND}_sorted_merged_rmdups.bam
BAI=${IND}_sorted_merged_rmdups.bam.bai


# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=samtools_env
# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate


echo "Processing chromosome ${CHR}"
samtools view -bh ${BAM} ${CHR} > ${IND}_${CHR}_sorted_merged_rmdups.bam
samtools index ${IND}_${CHR}_sorted_merged_rmdups.bam
