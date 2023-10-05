#! /bin/bash
#eval "$(conda shell.bash hook)"
#conda activate py37
module load picard
module load bcftools

# batch genotype from RNAseq SNP calls, imputed genotype data
cp /hpcdata/sg/sg_data/illumina_NCI_runs/PregCoV_NICA_PregCoV_bulk/bam/PregCoV/hg19_snp150_markdup_rmi.q4.bcf hg19_snp150_markdup_rmi.q4.bcf
bcftools index hg19_snp150_markdup_rmi.q4.bcf

# decompress
bcftools view hg19_snp150_markdup_rmi.q4.bcf -O v -o PregCoV.Demux.vcf

# get sample names
bcftools query -l PregCoV.Demux.vcf
bcftools query -l PregCoV.Demux.vcf > vcfsamples.txt

#edit the vcfsamples.txt to cleanup the sample names
grep -oE "[[:upper:]]+_PREG_[[:digit:]]+|[[:upper:]]+_STA_[[:digit:]]+|PREG_[[:digit:]]|CHI[[:digit:]]+c|CHI[[:digit:]]+|CHI[[:digit:]]+rc" vcfsamples.txt > vcfsamples.rename.txt
bcftools reheader PregCoV.Demux.vcf -s vcfsamples.rename.txt -o PregCoV.Demux.reheaded.vcf
bcftools query -l PregCoV.Demux.reheaded.vcf

bcftools stats -s- PregCoV.Demux.reheaded.vcf > bcfstats.output.vchk
bgzip PregCoV.Demux.reheaded.vcf
bcftools index PregCoV.Demux.reheaded.vcf.gz

#filtering
bcftools view PregCoV.Demux.reheaded.vcf.gz \
	--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > PregCoV.Demux.reheaded.sublex.vcf

# modify vcf header to work with the demuxlet pipeline
awk '{print $1}' PregCoV.Demux.reheaded.sublex.vcf | uniq
bcftools view -h PregCoV.Demux.reheaded.sublex.vcf > vcfheader
# manually edit the vcfheader
#bcftools reheader PregCoV.Demux.reheaded.sublex.vcf -h vcfheader.reorder -o PregCoV.Demux.CHI014.sublex.vcf
#awk '{print $1}' PregCoV.Demux.CHI014.sublex.vcf | uniq

# filter to sites with variants
#bcftools view -i 'AC>0' PregCoV.Demux.reheaded.sublex.vcf > PregCoV.fordemux.w.Variants.vcf
# add AF (Allele frequency) field
#bcftools +fill-tags PregCoV.fordemux.w.Variants.vcf -- -t AF | bcftools view -O z5 > PregCoV.fordemux.w.Variants.AF.vcf.gz

#filter to exclude ENCODE blacklisted regions for each batch
module load vcftools/0.1.17-jydu-goolf-1.7.20-Perl-5.22.2
vcftools --gzvcf PregCoV.Demux.reheaded.sublex.vcf --exclude-bed ENCFF001TDO.bed --recode --out PregCoV.fordemux
vcftools --vcf PregCoV.fordemux.recode.vcf --minGQ 10 --max-missing 1.0 --recode --out PregCoV.fordemux.filt
awk '{print $1}' PregCoV.fordemux.filt.recode.vcf | uniq
   
# for demuxlet
#bgzip PregCoV.fordemux.filt.recode.vcf
#bcftools index PregCoV.fordemux.filt.recode.vcf.gz
#bcftools view PregCoV.fordemux.filt.recode.vcf.gz \
#	--regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > PregCoV.fordemux.lexographic.vcf 

#awk '{print $1}' PregCoV.fordemux.lexographic.vcf | uniq
#bcftools query -l PregCoV.fordemux.lexographic.vcf

# change chr annotations to match bam files if necessary
bcftools annotate --rename-chrs rename_chrs.txt PregCoV.fordemux.filt.recode.vcf > PregCoV.fordemux.lexographic.wo.chr.prefix.vcf
grep '^#' PregCoV.fordemux.lexographic.wo.chr.prefix.vcf > PregCoV.fordemux.lexographic.wo.chr.prefix.ordered.vcf && grep -v '^#' PregCoV.fordemux.lexographic.wo.chr.prefix.vcf | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> PregCoV.fordemux.lexographic.wo.chr.prefix.ordered.vcf

# PCA
VCF="PregCoV.fordemux.lexographic.wo.chr.prefix.ordered.vcf"
/sysapps/cluster/software/PLINK/2.00-x86_64-alpha/plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:#[b37]\$r,\$a --indep-pairwise 50 10 0.1 --out plink --new-id-max-allele-len 31
/sysapps/cluster/software/PLINK/2.00-x86_64-alpha/plink2 --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:#[b37]\$r,\$a --extract plink.prune.in --make-bed --pca --out plink --new-id-max-allele-len 31
