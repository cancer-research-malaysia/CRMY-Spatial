setwd("~/ssgsea")

## load libraries
library(GSVA)

#GeneSet <- read.table("Example_gene_sets_GSVA.txt", header=TRUE) # Need to add NAs to make lists equal length for read.table
GeneSet <- read.table("Genesets_for_GSVA.txt", sep="\t", header=TRUE)
GeneSet <- as.list(GeneSet)
GeneSet <- lapply(GeneSet, function(x) x[!is.na(x)]) # remove NAs

#"SD507",
for (sample in c("SD560","SD1043", "SD1182","SD1225", "SD683","SD693","SD781")) {
  
ExprSet <- read.csv(paste("/stereoseq/all_samples/tumour_gene_counts/",sample,".tumour.gene_counts.csv", sep=""), 
                    sep=",", header=TRUE, row.names = 1)
ExprSet <- t(ExprSet)

ExprMat <- as.matrix(ExprSet[,-1])
# rownames(ExprMat)<-ExprSet[,1]


GSVA_out <- gsva(ExprMat, GeneSet, kcdf="Gaussian", verbose=T, parallel.sz=4)


GsVA_res_df <- t(GSVA_out)
GsVA_res_df <- as.data.frame(GsVA_res_df)
GsVA_res_df <- format(round(GsVA_res_df, 4), nsmall = 4)


GsVA_res_df <- cbind(row.names(GsVA_res_df),GsVA_res_df) 
colnames(GsVA_res_df)[1]<-c("bin")


#write.table(GsVA_res_df, file='Example_GSVA_out_test.txt', quote=FALSE, sep='\t', row.names = F)
write.table(GsVA_res_df, file=paste(sample,".tumour.gsva_out2.txt", sep=""), quote=FALSE, sep='\t', row.names = F)

}
