#seurat analysis

library(Seurat)
library(tidyverse)
library(Matrix)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

counts <- readMM(file = "../raw_data_full/seurat_analysis_data/master_matrix.txt")
genes <- read.csv("../raw_data_full/seurat_analysis_data/master_genes.csv")
metadata <- read.csv("../raw_data_full/metadata_final.csv")

metadata <- metadata[,4:8]

rownames(counts) <- genes$x
colnames(counts) <- metadata$cells

rownames(metadata) <- metadata$cells

#create original seurat object with counts data
#use the same QC as atlas study, if not mentioned I will assume default

seurat <- CreateSeuratObject(counts = counts, project = "male_vs_female_mice_spinal", meta.data = metadata)
seurat

#original seurat has 33992400 cells

seurat@meta.data$orig.ident <- as.factor(seurat@meta.data$study)
seurat@meta.data$study <- as.factor(seurat@meta.data$study)
seurat@meta.data$sex <- as.factor(seurat@meta.data$sex)
#metadata <- seurat@meta.data

#saveRDS(seurat, file = "../raw_data_full/my_data.rds")

#show percent of genes linked to mitochondria
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

#use the same QC as atlas study, if not mentioned I will assume default
#seurat <- subset(seurat, subset = nFeature_RNA > 200 & percent.mt < 5)
seurat <- subset(seurat, subset = nFeature_RNA > 100 & percent.mt < 5)

#after subset there's 244222 cells
#save for SeqSeek analysis.
saveRDS(seurat, file = "../raw_data_full/query_large.rds")

#NEXT STEP : Run through 3b, then Run Through SeqSeek pipeline
#this query object will fail the second step, this is expected behaviour.

#---------------------------------------------------------------------

#lognormalized to 10,000 counts.
#seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable features
#seurat <- FindVariableFeatures(seurat, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

#all.genes <- rownames(seurat)

#saveRDS(seurat, file = "../raw_data_full/my_data.rds")

#consider using SCTransform instead of scale data for this step, as recommended by authors
#need a lot of memory for this step

#seurat <- ScaleData(seurat, model.use = "linear", vars.to.regress = c("percent.mt", "nCount_RNA"), features = all.genes)

#start integration
#pre_split_seurat <- readRDS("../raw_data_full/pre_split_seurat.rds")

#split_seurat <- SplitObject(pre_split_seurat, split.by = "orig.ident")

#normalize using log normalization at a scale of 10 000
# find variable features with mean.var.plot selection method with custom cutoffs.
#split_seurat <- lapply(X = split_seurat, FUN = function(x) {
#  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
#  x <- FindVariableFeatures(x, selection.method = "mean.var.plot", mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5, Inf))
#})

#features <- SelectIntegrationFeatures(object.list = split_seurat)

#integrate using 20 PCs.
#seurat_anchors <- FindIntegrationAnchors(object.list = split_seurat, anchor.features = features, dims = 1:20)

# this command creates an 'integrated' data assay
#seurat_combined <- IntegrateData(anchorset = seurat_anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
#DefaultAssay(seurat_combined) <- "integrated"

#regress out percent.mt and the number of RNAs as covariates.
#seurat_combined_scale <- ScaleData(seurat_combined, model.use = "linear", vars.to.regress = c("percent.mt", "nCount_RNA") )

#find 100 PCS.
#PCA_final <- RunPCA(seurat_combined_scale, features = VariableFeatures(object = seurat_combined_scale), npcs = 100 )

#saveRDS(PCA_final, file = "pca_final.rds")

#seurat <- readRDS("pca_final.rds")

#ElbowPlot(seurat, ndims = 20)
#?ElbowPlot

#choose 15 dims.

#seurat <- FindNeighbors(seurat, dims = 1:15)
#seurat <- FindClusters(seurat, resolution = 0.8)

#seurat <- RunUMAP(seurat, dims = 1:15)

#markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
#?FindAllMarkers

#saveRDS(seurat, "data_done.rds")
#write.csv(markers, "cluster_markers_2.csv")

