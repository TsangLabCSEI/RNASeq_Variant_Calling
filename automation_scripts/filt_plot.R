library(vcfR)
library(stringr)
library(pheatmap)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

#read vcf
vcf <- read.vcfR(args[1], verbose = FALSE)
gt <- extract.gt(vcf, return.alleles = FALSE)
df_org <-sapply(as.data.frame(gt), function(x) {str_detect(x, "0/1") | str_detect(x, "1/0") | str_detect(x, "1/1")})
#convert NAs to 0
df_org[is.na(df_org)] <- 0
CHROMs <- getFIX(vcf)[, 'CHROM']

samples <- args[3:length(args)]
samples <- append(samples, "FORMAT", 0)
df <- df_org[,which(colnames(df_org) %in% samples)]

m <- as.data.frame(matrix(0, ncol = ncol(df), nrow = nrow(df)))
colnames(m) <- colnames(df)
onevar <- as.data.frame(matrix(0, ncol = ncol(df), nrow = nrow(df)))
for (i in 1:length(colnames(df))) {
  onevar[which((df[,i] == 1) & (rowSums(df)==1)), i] = 1
  idx = sample(which(onevar[,i] == 1), as.numeric(args[2]))
  m[idx,i] = 1
}

write.vcf(vcf[which(rowSums(m)==1), c(1,which(colnames(df_org) %in% samples)+1)],"filtered_variants.vcf.gz")

#free up space and load in subselected vcf to generate plots
vcf <- read.vcfR("filtered_variants.vcf.gz", verbose = FALSE)
gt <- extract.gt(vcf, return.alleles = FALSE)
df<-sapply(as.data.frame(gt), function(x) {str_detect(x, "0/1") | str_detect(x, "1/0") | str_detect(x, "1/1")})
#convert NAs to 0
df[is.na(df)] <- 0
CHROMs <- getFIX(vcf)[, 'CHROM']

#plot all variants for each patient
pheatmap(df,cluster_cols=FALSE,show_rownames = FALSE,show_colnames = FALSE, filename="./allvars.png")

#plot number of unique variants per patient
png("./numuniquevars_filtered.png", width = 800, height = 480)
par(mar=c(5,8,4,1)+.1)
barplot(colSums(df),las=2,horiz = TRUE)
dev.off()

#plot number of unique variants per chromosome per patient
df2 <- as.data.frame(df)
df2$CHROM <- CHROMs

#initialize matrices
c <- data.frame(matrix(0, length(colnames(df2))-1, 24))
colnames(c) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X", "Y", "MT")
rownames(c) <- colnames(df2[1:length(colnames(df2))-1])

#Count var per chromosome
for (i in 1:length(colnames(c))) {
  datset <- subset(df2, CHROM == colnames(c)[i])
  for (j in 1:(nrow(c))) {
    c[j,i] <- sum(datset[,j])
  }
}

#plot
pheatmap(as.matrix(c),cluster_cols=FALSE, filename="./uniqevarsbychromosome.png")
