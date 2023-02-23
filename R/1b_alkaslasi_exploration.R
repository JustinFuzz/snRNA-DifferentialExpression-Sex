#alkaslasi has 3 x 2 different matrices. 2 being male and female.
#we will atleast want to merge the matrices belonging to the same groups here, if possible.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Matrix)
library(dplyr)

cerv_f_dir = "../raw_data_full/alkaslasi/Cervical_F_raw_feature_bc_matrix/"
cerv_m_dir = "../raw_data_full/alkaslasi/Cervical_M_raw_feature_bc_matrix/"

lumb_f_dir = "../raw_data_full/alkaslasi/Lumbar_F_raw_feature_bc_matrix/"
lumb_m_dir = "../raw_data_full/alkaslasi/Lumbar_M_raw_feature_bc_matrix/"

thor_f_dir = "../raw_data_full/alkaslasi/Thoracic_F_raw_feature_bc_matrix/"
thor_m_dir = "../raw_data_full/alkaslasi/Thoracic_M_raw_feature_bc_matrix/"


load_matrix_F <- function(matrix_dir){
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = paste0(barcode.names$V1, "_F")
  rownames(mat) = feature.names$V1
  return(mat)
}

load_matrix_M <- function(matrix_dir){
  barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
  features.path <- paste0(matrix_dir, "features.tsv.gz")
  matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = paste0(barcode.names$V1, "_M")
  rownames(mat) = feature.names$V1
  return(mat)
}

cerv_f_mat <- load_matrix_F(cerv_f_dir)
cerv_m_mat <- load_matrix_M(cerv_m_dir)

lumb_f_mat <- load_matrix_F(lumb_f_dir)
lumb_m_mat <- load_matrix_M(lumb_m_dir)

thor_f_mat <- load_matrix_F(thor_f_dir)
thor_m_mat <- load_matrix_M(thor_m_dir)

#get the mean of each group for each sex.

mean_f_mat <- (thor_f_mat + lumb_f_mat + cerv_f_mat) / 3
mean_m_mat <- (thor_m_mat + lumb_m_mat + cerv_m_mat) / 3

#add unique cell identifier to each group.
add_identifer <- function(matrix, identifer){
  
  colnames <- colnames(matrix)
  colnames <- paste0(identifer, colnames)
  colnames(matrix) <- colnames
  
  return(matrix)
}

#cerv_f_mat <- add_identifer(cerv_f_mat, "cerv_f_")
#cerv_m_mat <- add_identifer(cerv_m_mat, "cerv_m_")

#lumb_f_mat <- add_identifer(lumb_f_mat, "lumb_f_")
#lumb_m_mat <- add_identifer(lumb_m_mat, "lumb_m_")

#thor_f_mat <- add_identifer(thor_f_mat, "thor_f_")
#thor_m_mat <- add_identifer(thor_m_mat, "thor_m_")

mean_f_mat <- add_identifer(mean_f_mat, "mean_f_")
mean_m_mat <- add_identifer(mean_m_mat, "mean_m_")


#column bind matrices that belong to either male or female groups, without modification.

#f_matrix <- cbind(thor_f_mat, cbind(cerv_f_mat,lumb_f_mat))
#m_matrix <- cbind(thor_m_mat, cbind(cerv_m_mat,lumb_m_mat))

full_matrix <- cbind(mean_f_mat, mean_m_mat)

#blum $ alkasli use ensemble annotations of all genes.... so we need to convert them.
#using biomart package from bioconductors

library(biomaRt)

#we are using mmusculus_gene_ensembl
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

#get master list of ensembl codes from alkasli
gene_list_alkas <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"), values = rownames(full_matrix), mart = ensembl)

forbidden_index <- which(gene_list_alkas$external_gene_name == "")

subset_matrix <- full_matrix[rownames(full_matrix) %in% gene_list_alkas$ensembl_gene_id,]

#keep order matching with the current subset matrix.
gene_list_alkas <- gene_list_alkas[match(rownames(subset_matrix), gene_list_alkas$ensembl_gene_id),]

#remove rows with empty gene names (from forbidden index)
subset_matrix <- subset_matrix[ !(1:nrow(subset_matrix) %in% forbidden_index),]

subset_gene_list <- gene_list_alkas[ !(rownames(gene_list_alkas) %in% forbidden_index),]

rownames(subset_matrix) <- subset_gene_list$external_gene_name
colnames(subset_matrix) <- colnames(full_matrix)

#save
write.csv(colnames(subset_matrix), file = "../raw_data_full/alkaslasi/mean_cells.csv")
write.csv(rownames(subset_matrix), file = "../raw_data_full/alkaslasi/mean_genes.csv")
writeMM(subset_matrix, file = "../raw_data_full/alkaslasi/mean_matrix_alkaslasi.txt")
