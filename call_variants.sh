#! /bin/bash

ml load SAMtools/1.11-GCCcore-10.2.0
ml load BCFtools/1.16-GCCcore-10.2.0

for i in {1..22};do
	samtools mpileup -A -q 4 -t AD,DP -l ../BED/chr_${i}.bed -f ../nfcore_rnavar/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa *.bam -g | bcftools call -m - -O v -f GQ >  allsnps_chr${i}.vcf;
done
