#source("https://bioconductor.org/biocLite.R")
#biocLite("GSVA")
setwd("~/ssgsea")

## load libraries
library(GSVA)


#GeneSet <- read.table("Example_gene_sets_GSVA.txt", header=TRUE) # Need to add NAs to make lists equal length for read.table
GeneSet <- read.table("gsva_pathways.txt", sep="\t", header=TRUE)
GeneSet <- as.list(GeneSet)
GeneSet <- lapply(GeneSet, function(x) x[!is.na(x)]) # remove NAs



#ExprSet <- read.table("Example_GSVA_voomed.txt", sep="\t", header=TRUE)
ExprSet <- read.table("/stereoseq/all_samples/rnaseq/MyBrCa_Batch1_and_2_counts_voomed.txt", sep="\t", header=TRUE)


ExprMat <- as.matrix(ExprSet[,-1])
rownames(ExprMat)<-ExprSet[,1]



GSVA_out <- gsva(ExprMat, GeneSet, kcdf="Gaussian", verbose=T, parallel.sz=4)


GsVA_res_df <- t(GSVA_out)
GsVA_res_df <- as.data.frame(GsVA_res_df)
GsVA_res_df <- format(round(GsVA_res_df, 4), nsmall = 4)

GsVA_res_df <- cbind(row.names(GsVA_res_df),GsVA_res_df) 
colnames(GsVA_res_df)[1]<-c("ID")

#write.table(GsVA_res_df, file='Example_GSVA_out_test.txt', quote=FALSE, sep='\t', row.names = F)
write.table(GsVA_res_df, file='MyBrCa_TNBC.gsva_out.txt', quote=FALSE, sep='\t', row.names = F)

