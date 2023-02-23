library(Seurat)
library(tidyverse)
library(Matrix.utils)
library(dplyr)
library(Matrix)
library(purrr)
library(tibble)
library(SingleCellExperiment)
library(apeglm)
library(DESeq2)
library(DEGreport)
library(magrittr)
library(sva)

#to be ran after updated SeqSeek Pipeline. With updated 'query.rds' object.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

seurat_small <- readRDS("results_small/query.rds")
seurat_large <- readRDS("results_large/query.rds")

seurat_small <- SetIdent(seurat_small, value = seurat_small@meta.data$predicted.id)

seurat_small <- subset(seurat_small, idents = c("Junk", "Doublets"), invert = TRUE)

seurat_large@meta.data$predicted.id <- query_small@meta.data$predicted.id

seurat <- seurat_large
rm(seurat_large)
rm(seurat_small)

#following https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

#get raw counts
counts <- seurat@assays$RNA@counts
metadata <- seurat@meta.data
metadata$sample_id <- metadata$sample

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- as.factor(seurat@active.ident)
metadata$study <- as.factor(metadata$study)
metadata$group_id <- metadata$

levels(metadata_clean$study)
levels(metadata_clean$study) <- c("alkaslasi.blum", "alkaslasi.blum", "sathyamurthy")

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

colData(sce)
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))

# Total number of samples 
ns <- length(sids)
ns

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei

#QC was already done...

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)

pb[1:6, 1:6]

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))

get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()

# Get cluster IDs for each of the samples

samples_list <- map(1:length(kids), get_sample_ids)

get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

gg_df <- data.frame(cluster_id = de_cluster_ids,
                    sample_id = de_samples)

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id", "study")]) 

gg_df$sample_id <- as.factor(gg_df$sample_id)

metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id, study) 

metadata   

# Generate vector of cluster IDs
metadata$cluster_id[metadata$cluster_id == "Peripheral Glia/BC Derivative"] <- "Peripheral Glia.BC Derivative"

clusters <- levels(as.factor(metadata$cluster_id))
clusters

# Likelihood ratio test
dir.create("DESeq2/results")

# DESeq2
get_dds_LRTresults <- function(x){
  
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  counts <- pb[[clusters[x]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  
  rownames(cluster_metadata) <- colnames(cluster_counts)
  
  print(all(rownames(cluster_metadata) == colnames(cluster_counts))  )      
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ study + group_id)
  
  #use sva to reduce unnaccounted for variation.
  dat  <- counts(dds, normalized = FALSE)
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  mod  <- model.matrix(~ study + group_id, colData(dds))
  mod0 <- model.matrix(~   1, colData(dds))
  
  svseq <- svaseq(dat, mod, mod0)
  
  ddssva <- dds
  design <- ~ study + group_id
  
  for (i in ncol(svseq$sv):1){
    number <- paste0("SV",i)
    ddssva$number <- svseq$sv[,i]
    update(design, ~ number + .)
  }

  design(ddssva) <- design
  
  dds_lrt <- DESeq(ddssva, test="Wald")
  
  # Extract results
  res_LRT <- results(dds_lrt)
  
  # Create a tibble for LRT results
  res_LRT_tb <- res_LRT %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
  # Save all results
  write.csv(res_LRT_tb,
            paste0("DESeq2/wald/", clusters[x], "_wald_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Subset to return genes with padj < 0.05
  sigLRT_genes <- res_LRT_tb %>% 
    filter(padj < 0.05)
  
  # Save sig results
  write.csv(sigLRT_genes,
            paste0("DESeq2/wald/", clusters[x], "_wald_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
}

map(1:length(clusters), get_dds_LRTresults)

#now run a custom DESeq2 analysis. 1. all cells 2. SDH neurons. 3. DDH neurons. 4. SDH + DDH neurons
#do this with fresh metadata
doublets <- metadata$cells[metadata$cell_type == "Doublets"]
junk <- metadata$cells[metadata$cell_type == "Junk"]

bad_cells <- c(doublets, junk)

metadata <- gg_df %>%
  dplyr::select(cluster_id, sample_id, group_id, study)

metadata$cluster_id <- as.factor(metadata$cluster_id)

metadata_clean <- metadata[which(!metadata$cluster_id == "Junk"),]
metadata_clean <- metadata_clean[which(!metadata_clean$cluster_id == "Doublets"),]

droplevels(metadata_clean$cluster_id)
#metadata_clean <- data.frame(sample_id = metadata_clean$sample_id, group_id = metadata_clean$group_id, study_id = metadata_clean$study)

counts_clean <- counts[, !colnames(counts) %in% bad_cells]
#metadata_clean <- metadata[!metadata$cluster_id %in% c("Doublets", "Junk"),]


#make sce
sce_clean <- SingleCellExperiment(assays = list(counts = counts_clean), 
                            colData = metadata_clean)

groups_clean <- colData(sce_clean)[, c("sample_id")]

# Aggregate across cluster-sample groups
pb_clean <- aggregate.Matrix(t(counts(sce_clean)), 
                       groupings = groups_clean, fun = "sum")

pb_real <- t(pb_clean)

metadata_clean <- ei

dds_clean <- DESeqDataSetFromMatrix(pb_real,
                              colData = metadata_clean, 
                              design = ~ study + group_id)

dat  <- counts(dds_clean, normalized = FALSE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ study + group_id, colData(dds_clean))
mod0 <- model.matrix(~   1, colData(dds_clean))

svseq <- svaseq(dat, mod, mod0)

ddssva <- dds_clean
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
design(ddssva) <- ~ SV1 + SV2 + SV3 + study + group_id

ddssva_wald <- DESeq(ddssva, test="Wald")

rldsva <- rlog(ddssva_wald, blind=TRUE)

res <- results(ddssva_wald)

# Create a tibble for LRT results
res <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Save all results
write.csv(res,
          paste0("cells_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

# Subset to return genes with padj < 0.05
sigLRT_genes <- res %>% 
  filter(padj < 0.05)

# Save sig results
write.csv(sigLRT_genes,
          paste0("cells_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


# 2. 

#ddh neurons = 


