#remap genome if necessary...

library(Seurat)
library(dplyr)
library(readr)
library(Matrix)
library(MatrixExtra)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#need to find genome of all 3 experiments
#need to make sure we're using the same gene list as them (25419 genes)

#sathyamurthy: Genome_build: mm10 (GRCm38)
#blum: Genome_build: mm10 (mm10-3.0.0-premrna)
#alkaslasi: Genome_build: mm10

#-----------------------------------------------------------------
#new counts

counts <- readMM(file = "../raw_data_full/master_matrix.txt")
genes <- read.csv("../raw_data_full/master_genes.csv")
metadata <- read.csv("../raw_data_full/metadata_final.csv")

#need to downlaod this file from the SeqSeek github under 'model'
current_genes <- read.csv(file='../raw_data_full/filter_NN_genes.tsv', header = FALSE)

rownames(counts) <- genes$x
colnames(counts) <- metadata$cells

reorder_idx <- match(current_genes$V1,rownames(counts))
reorder_idx[is.na(reorder_idx)] <- 0 #replace with 0 because sparsematrix doesn't like nas

counts <- counts[reorder_idx,]

which(!current_genes$V1 %in% rownames(counts))

add_missing_genes <- function(smat, master_list){
  
  #create a sparse matrix of empty genes that are missing, then rbind to current matrix.
  missing_genes <- setdiff(master_list, rownames(smat))
  
  #create matrix with all 0's and make it a sparse matrix.
  
  #need to create a sparse matrix directly here instead of a matrix.
  missing_smat <- emptySparse(nrow = length(missing_genes), ncol = length(colnames(smat)), format = "C", dtype = "d" )
  
  colnames(missing_smat) <- colnames(smat)
  rownames(missing_smat) <- missing_genes
  
  #bind the two matrices together
  full_smat <- rbind(smat, missing_smat)
  
  #change order to match master list
  full_smat <- full_smat[master_list,]
  
  return(full_smat)
}

#missing 1101 genes, add them to the bottom of the matrix with 0 then??
counts <- add_missing_genes(counts, current_genes$V1)

rownames(metadata) <- metadata$cells

#create second seurat object
seurat2 <- CreateSeuratObject(counts = counts, project = "male_vs_female_mice_spinal", meta.data = metadata)

seurat2@meta.data$orig.ident <- as.factor(seurat2@meta.data$study)
seurat2@meta.data$study <- as.factor(seurat2@meta.data$study)
seurat2@meta.data$sex <- as.factor(seurat2@meta.data$sex)

#show percent of genes linked to mitochondria
seurat2[["percent.mt"]] <- PercentageFeatureSet(seurat2, pattern = "^MT-")

#use the same QC as atlas study, if not mentioned I will assume default
seurat2 <- subset(seurat2, subset = nFeature_RNA > 100 & percent.mt < 5)
saveRDS(seurat2, file = "../raw_data_full/query_small.rds")

#now you can run SeqSeek with both query seurat objects.
# you just need to change their name to 'query.rds' and put them in data.
# but just follow their github: https://github.com/ArielLevineLabNINDS/SeqSeek_Classify_Full_Pipeline