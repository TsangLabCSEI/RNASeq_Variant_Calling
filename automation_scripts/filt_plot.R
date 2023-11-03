library(vcfR)
library(stringr)
library(pheatmap)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

#read vcf
vcf <- read.vcfR(args[1], verbose = FALSE)
gt <- extract.gt(vcf, return.alleles = FALSE)
df<-sapply(as.data.frame(gt), function(x) {str_detect(x, "0/1") | str_detect(x, "1/0") | str_detect(x, "1/1")})
#convert NAs to 0
df[is.na(df)] <- 0
CHROMs <- getFIX(vcf)[, 'CHROM']

#write vcf with unique variants for each patient
#write.vcf(vcf[which(rowSums(df)==1),],"filtered_variants.vcf.gz")

#plot all variants for each patient
pheatmap(df,cluster_cols=FALSE,show_rownames = FALSE,show_colnames = FALSE, filename="./allvars.png")

#plot unique variants for each patient
pheatmap(df[which(rowSums(df)==1),],cluster_cols=FALSE,show_rownames = FALSE,show_colnames = TRUE, filename="./uniquevars.png")

#plot number of unique variants per patient
png("./numuniquevars.png", width = 800, height = 480)
par(mar=c(5,8,4,1)+.1)
barplot(colSums(df[which(rowSums(df)==1),]),las=2,horiz = TRUE)
dev.off()

#plot number of unique variants per chromosome per patient
df2 <- as.data.frame(df)
df2$CHROM <- CHROMs
df2 <- df2[which(rowSums(df)==1), ]

c <- data.frame(matrix(0, length(colnames(df)), 24))
colnames(c) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X", "Y", "MT")
rownames(c) <- colnames(df)

for (i in 1:length(colnames(c))) {
  datset <- subset(df2, CHROM == toString(i))
  for (j in 1:length(rownames(c))) {
    c[j,i] <- sum(datset[j])
  }
}

pheatmap(as.matrix(c),cluster_cols=FALSE, filename="./uniqevarsbychromosome.png")
