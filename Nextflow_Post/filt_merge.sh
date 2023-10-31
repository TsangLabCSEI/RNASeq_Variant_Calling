#!/bin/sh

bcftools view -i $2 $1 > ${1%.haplotypecaller.filtered.vcf}.DP2000.vcf

#make sure bed file is same folder or the path is correct       
vcftools --vcf ${1%.haplotypecaller.filtered.vcf}.DP2000.vcf --exclude-bed ENCFF001TDO.bed --recode --out ${1%.haplotypecaller.filtered.vcf}.DP2000.fordemux
vcftools --vcf ${1%.haplotypecaller.filtered.vcf}.DP2000.fordemux.recode.vcf --minGQ 10 --max-missing 1.0 --recode --out ${1%.haplotypecaller.filtered.vcf}.DP2000.fordemux.filt

bgzip ${1%.haplotypecaller.filtered.vcf}.DP2000.fordemux.filt.recode.vcf
bcftools index ${1%.haplotypecaller.filtered.vcf}.DP2000.fordemux.filt.recode.vcf.gz
