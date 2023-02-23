#BiocManager::install("rhdf5")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rhdf5)
library(Matrix)
library(tibble)

female_5_6 <- h5read("../raw_data_full/blum/GSM4911292_female5_6_25.h5", "matrix")
female_5_7 <- h5read("../raw_data_full/blum/GSM4911294_female5_7_30.h5", "matrix")
male_4_4 <- h5read("../raw_data_full/blum/GSM4911293_male4_pilot_4_25.h5", "matrix")

createSparseMatrix_F <- function(h5){
  
  my_barcodes <- as.character(h5[["barcodes"]])
  
  counts <- sparseMatrix(
    dims = h5$shape,
    i = as.numeric(h5$indices),
    p = as.numeric(h5$indptr),
    x = as.numeric(h5$data),
    index1 = FALSE
  )
  colnames(counts) <- paste0(my_barcodes, "_F")
  rownames(counts) <- as.data.frame(h5[["features"]])$id
  return(counts)
}

createSparseMatrix_M <- function(h5){
  
  my_barcodes <- as.character(h5[["barcodes"]])
  
  counts <- sparseMatrix(
    dims = h5$shape,
    i = as.numeric(h5$indices),
    p = as.numeric(h5$indptr),
    x = as.numeric(h5$data),
    index1 = FALSE
  )
  colnames(counts) <- paste0(my_barcodes, "_M")
  rownames(counts) <- as.data.frame(h5[["features"]])$id
  return(counts)
}

matrix_female_5_6 <- createSparseMatrix_F(female_5_6)
matrix_female_5_7 <- createSparseMatrix_F(female_5_7)
matrix_male_4_4 <- createSparseMatrix_M(male_4_4)

#make unique cells names for each.
colnames <- colnames(matrix_female_5_6)
colnames <- paste0("5_6_", colnames)
colnames(matrix_female_5_6) <- colnames

colnames <- colnames(matrix_female_5_7)
colnames <- paste0("5_7_", colnames)
colnames(matrix_female_5_7) <- colnames

colnames <- colnames(matrix_male_4_4)
colnames <- paste0("4_4_", colnames)
colnames(matrix_male_4_4) <- colnames

#full_matrix <- cbind(matrix_male_4_4, cbind(matrix_female_5_6, matrix_female_5_7))

#blum $ alkasli use ensemble annotations of all genes.... so we need to convert them.
#using biomart package from bioconductor

library(biomaRt)

listDatasets(useMart("ensembl"))
listFilters(ensembl)

#we are using mmusculus_gene_ensembl
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", mirror = "useast")

#get master list of ensembl codes from and blum

convert_gene_names <- function(full_matrix){

  #get list of correct genes needed
  gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "external_gene_name"), values = rownames(full_matrix), mart = ensembl)
  
  #get index of empty genes
  forbidden_index <- which(gene_list$external_gene_name == "")
  
  #subset matrix that contains the genes that were associated with the search
  subset_matrix <- full_matrix[rownames(full_matrix) %in% gene_list$ensembl_gene_id,]
  
  #match the order of gene names to the matrix
  gene_list <- gene_list[match(rownames(subset_matrix), gene_list$ensembl_gene_id),]
  
  #remove rows with empty gene names (from forbidden index)
  subset_matrix <- subset_matrix[ !(1:nrow(subset_matrix) %in% forbidden_index),]
  
  subset_gene_list <- gene_list[ !(rownames(gene_list) %in% forbidden_index),]
  
  rownames(subset_matrix) <- subset_gene_list$external_gene_name
  colnames(subset_matrix) <- colnames(full_matrix)
  
  return(subset_matrix)
  
}

matrix_female_5_6 <- convert_gene_names(matrix_female_5_6)
matrix_female_5_7 <- convert_gene_names(matrix_female_5_7)
matrix_male_4_4 <- convert_gene_names(matrix_male_4_4)

#save
save_information <- function(matrix, name){
  new_dir <- paste0("../raw_data_full/blum/", paste0(name, "/") )
  dir.create(new_dir)
  
  write.csv(colnames(matrix), file = paste0( new_dir, "cells.csv") )
  write.csv(rownames(matrix), file = paste0( new_dir, "genes.csv") )
  writeMM(matrix, file = paste0(new_dir, paste0(name, ".txt")) )
}

#save_information(matrix_female_5_6, "matrix_female_5_6")
#save_information(matrix_female_5_7, "matrix_female_5_7")
#save_information(matrix_male_4_4, "matrix_male_4_4")

subset_matrix <- cbind(matrix_female_5_6, cbind(matrix_female_5_7, matrix_male_4_4))

write.csv(colnames(subset_matrix), file = "../raw_data_full/blum/cells.csv")
write.csv(rownames(subset_matrix), file = "../raw_data_full/blum/genes.csv")
writeMM(subset_matrix, file = "../raw_data_full/blum/matrix_blum.txt")


