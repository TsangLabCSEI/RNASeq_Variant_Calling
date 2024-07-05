#!/bin/sh


bcftools view -i $2 $1 > ${1%.vcf}.DP${3}.vcf

#make sure bed file is same folder or the path is correct       
vcftools --vcf ${1%.vcf}.DP${3}.vcf --exclude-bed ENCFF001TDO.bed --recode --out ${1%.vcf}.DP${3}.fordemux
vcftools --vcf ${1%.vcf}.DP${3}.fordemux.recode.vcf --minGQ 10 --max-missing 1.0 --recode --out ${1%.vcf}.DP${3}.fordemux.filt

bgzip ${1%.vcf}.DP${3}.fordemux.filt.recode.vcf
bcftools index ${1%.vcf}.DP${3}.fordemux.filt.recode.vcf.gz
