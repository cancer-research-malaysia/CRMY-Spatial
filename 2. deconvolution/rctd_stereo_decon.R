########################## RCTD deconvolution for Stereo-seq data ##########################
## load libraries
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(spacexr)
#library(scrattch.io)
library(Matrix)
#library(spam)
#library(spam64)

########################## create reference ##########################
## reference from: https://doi.org/10.1038/s41588-021-00911-1 

ref_name = "GSE176078"
expression_matrix <- ReadMtx(
  mtx = paste("reference/scRNAseq/",ref_name,"/matrix.mtx.gz", sep=""),
  features = paste("reference/scRNAseq/",ref_name,"/features.tsv.gz", sep=""),
  cells = paste("reference/scRNAseq/",ref_name,"/barcodes.tsv.gz", sep="")
)
reference_matrix <- CreateSeuratObject(counts = expression_matrix)

ref_counts <- reference_matrix@assays[["RNA"]]@layers[["counts"]] #cells@assays[["RNA"]]@counts
rownames(ref_counts) <- rownames(reference_matrix@assays[["RNA"]])
colnames(ref_counts) <- colnames(reference_matrix@assays[["RNA"]])

for (ref_sample in list.files("reference/scRNAseq/GSE176078/tnbc/")) {
  df1 <- read.csv(paste("reference/scRNAseq/GSE176078/tnbc/",ref_sample,"/metadata.csv", sep=""),
                  header=TRUE, row.names="X")
  if (exists("dfs")) {
    dfs <- rbind(dfs, df1)
  }
  else {
    dfs <- df1
  }
}
celltype_subset <- dfs[,"celltype_subset"]
# celltype_minor <- dfs[,"celltype_minor"]
df <- as.data.frame(celltype_subset)
# df <- as.data.frame(celltype_minor)
cell_types <- factor(as.matrix(df))
names(cell_types) <- rownames(dfs)

reference <- Reference(ref_counts, cell_types)

########################## create spatialRNA ##########################

bin_size <- 30
#sample_name <- "SD507"
for (sample_name in c("SD507","SD560","SD1043", # A3B del/del
                      "SD1182","SD1225", # A3B del/WT
                      "SD683","SD693","SD781")) {
  
  #reference <- Reference(ref_counts, cell_types)
  data_path <- paste("data/bin", bin_size, "/", 
                     sample_name, ".bin", bin_size, ".seurat", sep="")
  if (!file.exists(paste(data_path, ".h5seurat", sep=""))) {
    print("Converting H5Seurat")
    Convert(paste(data_path, ".h5ad", sep=""), #_processed
            dest="h5seurat", assay="RNA", overwrite=TRUE)
  }
  print("Loading H5Seurat")
  adata <- LoadH5Seurat(paste(data_path, ".h5seurat", sep=""), #seurat_processed
                        assays="RNA", meta.data=FALSE, misc=FALSE)
  adata
  
  coords <- as.data.frame(adata@reductions[["spatial"]]@cell.embeddings) #as.data.frame(
  counts <- adata@assays[["RNA"]]@counts
  #spamcounts <- as.dgCMatrix.spam(as.spam(counts)) # use spam to handle large matrix
  #sparsecounts <- Matrix(counts, sparse=T)
  rawcounts <- t(readMM(paste("data/bin", bin_size, "/", 
                              sample_name,".bin",bin_size,".raw.counts.int.mtx", sep="")))
  rownames(rawcounts) <- rownames(counts)
  colnames(rawcounts) <- colnames(counts)
  spaceRNA <- SpatialRNA(coords, rawcounts) #, require_int=FALSE)
  
  # spaceRNA <- read.SpatialRNA(paste("data/bin", bin_size, sep=""), 
  #                             count_file = adata@assays[["RNA"]]@counts,
  #                             #count_file = paste(sample_name, ".bin", bin_size, ".processed.counts.csv", sep=""),
  #                             coords_file = adata@reductions[["spatial"]]@cell.embeddings
  #                             #coords_file = paste(sample_name, ".bin", bin_size, ".positions.csv", sep="")
  #                             #nUMI=NULL, require_int=FALSE
  #                             )
  
  ## clear memory
  rm(list=setdiff(ls(), c("spaceRNA","bin_size","reference","sample_name","ref_name")))
  
  ## run RCTD
  RCTD <- create.RCTD(spaceRNA, reference, CELL_MIN_INSTANCE=19)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  
  write.csv(RCTD@results, paste(sample_name,".bin", bin_size, ".RCTD.",ref_name,".full.results2.csv", sep=""))
}
