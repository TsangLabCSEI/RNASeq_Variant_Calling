library(vcfR)
library(stringr)
library(pheatmap)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

#read vcf
vcf <- read.vcfR(args[1], verbose = FALSE)
gt <- extract.gt(vcf, return.alleles = FALSE)
df<-sapply(as.data.frame(gt), function(x) {str_detect(x, "0/1") | str_detect(x, "1/0") | str_detect(x, "1/1")})
#convert NAs to 0
df[is.na(df)] <- 0
CHROMs <- getFIX(vcf)[, 'CHROM']

samples <- args[3:length(args)]

m <- as.data.frame(matrix(0, ncol = ncol(df), nrow = nrow(df)))
colnames(m) <- colnames(df)
onevar <- as.data.frame(matrix(0, ncol = ncol(df), nrow = nrow(df)))
for (i in 1:length(colnames(df))) {
  onevar[which((df[,i] == 1) & (rowSums(df)==1)), i] = 1
  idx = sample(which(onevar[,i] == 1), as.numeric(args[2]))
  m[idx,i] = 1
}

write.vcf(vcf[which(rowSums(m)==1), samples],"filtered_variants.vcf.gz")

#plot all variants for each patient
pheatmap(df,cluster_cols=FALSE,show_rownames = FALSE,show_colnames = FALSE, filename="./allvars.png")

#plot unique variants for each patient
pheatmap(df[which(rowSums(df)==1), samples],cluster_cols=FALSE,show_rownames = FALSE,show_colnames = TRUE, filename="./uniquevars.png")
pheatmap(df[which(rowSums(m)==1), samples],cluster_cols=FALSE,show_rownames = FALSE,show_colnames = TRUE, filename="./uniquevars_filtered.png")

#plot number of unique variants per patient
png("./numuniquevars.png", width = 800, height = 480)
par(mar=c(5,8,4,1)+.1)
barplot(colSums(df[which(rowSums(df)==1), samples]),las=2,horiz = TRUE)
dev.off()

png("./numuniquevars_filtered.png", width = 800, height = 480)
par(mar=c(5,8,4,1)+.1)
barplot(colSums(df[which(rowSums(m)==1), samples]),las=2,horiz = TRUE)
dev.off()

#plot number of unique variants per chromosome per patient
df2 <- as.data.frame(df)
df2$CHROM <- CHROMs
df2_filt <- as.data.frame(df2)[which(rowSums(m[,samples])==1), samples]
df2_filt$CHROM <- CHROMs[which(rowSums(m[,samples])==1)]
df2 <- df2[(which(rowSums(df)==1)),]

#initialize matrices
c <- data.frame(matrix(0, length(colnames(df2))-1, 24))
colnames(c) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X", "Y", "MT")
rownames(c) <- colnames(df2[1:length(colnames(df2))-1])

c_filt <- data.frame(matrix(0, length(colnames(df2_filt))-1, 24))
colnames(c_filt) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X", "Y", "MT")
rownames(c_filt) <- colnames(df2_filt[1:length(colnames(df2_filt))-1])

#Count var per chromosome
for (i in 1:length(colnames(c))) {
  datset <- subset(df2, CHROM == colnames(c)[i])
  for (j in 1:(nrow(c))) {
    c[j,i] <- sum(datset[,j])
  }
}

for (i in 1:length(colnames(c_filt))) {
  datset <- subset(df2_filt, CHROM == colnames(c_filt)[i])
  for (j in 1:nrow(c_filt)) {
    c_filt[j,i] <- sum(datset[,j])
  }
}

#plot
pheatmap(as.matrix(c),cluster_cols=FALSE, filename="./uniqevarsbychromosome.png")
pheatmap(as.matrix(c_filt),cluster_cols=FALSE, filename="./uniqevarsbychromosome_filtered.png")
