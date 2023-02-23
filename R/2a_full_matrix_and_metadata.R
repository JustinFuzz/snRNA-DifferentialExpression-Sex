#load 3 matrices.

library(Matrix)
library(MatrixExtra)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#they don't have dim names, so we need to attach the cells and genes in proper order for each
#before continuing.
sathy <- readMM("../raw_data_full/sathyamurthy/matrix_sathyamurthy.txt")
#haring <- readMM("../raw_data_full/haring/matrix_haring.txt")
blum <- readMM("../raw_data_full/blum/matrix_blum.txt")
alkas <- readMM("../raw_data_full/alkaslasi/mean_matrix_alkaslasi.txt")

#load cells and genes csv'

sathy_cells <- read.csv("../raw_data_full/sathyamurthy/cells.csv")
sathy_genes <- read.csv("../raw_data_full/sathyamurthy/genes.csv")

#haring_cells <- read.csv("../raw_data_full/haring/cells.csv")
#haring_genes <- read.csv("../raw_data_full/haring/genes.csv")

blum_cells <- read.csv("../raw_data_full/blum/cells.csv")
blum_genes <- read.csv("../raw_data_full/blum/genes.csv")

alkas_cells <- read.csv("../raw_data_full/alkaslasi/mean_cells.csv")
alkas_genes <- read.csv("../raw_data_full/alkaslasi/mean_genes.csv")

#add info to proper matrix.

colnames(sathy) <- sathy_cells$x
rownames(sathy) <- sathy_genes$x

#colnames(haring) <- haring_cells$x
#rownames(haring) <- haring_genes$x

colnames(blum) <- blum_cells$x
rownames(blum) <- blum_genes$x

colnames(alkas) <- alkas_cells$x
rownames(alkas) <- alkas_genes$x

#master_gene_list <- unique(c(sathy_genes$x, haring_genes$x, blum_genes$x, alkas_genes$x))
master_gene_list <- unique(c(sathy_genes$x, blum_genes$x, alkas_genes$x))


#which(alkas_genes$x == "")

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

sathy_full_mat <- add_missing_genes(sathy, master_gene_list)
#haring_full_mat <- add_missing_genes(haring, master_gene_list)
blum_full_mat <- add_missing_genes(blum, master_gene_list)
alkas_full_mat <- add_missing_genes(alkas, master_gene_list)

#master_matrix <- cbind(sathy_full_mat, cbind(haring_full_mat, cbind(blum_full_mat, alkas_full_mat) ) )
master_matrix <- cbind(sathy_full_mat, cbind(blum_full_mat, alkas_full_mat) )

#save 'master' matrix
writeMM(master_matrix, file = "../raw_data_full/master_matrix.txt")

#save genes as well.
write.csv(rownames(master_matrix), file = "../raw_data_full/master_genes.csv")

#Now, start creating metadata, then concatenate all matrices together into one large master matrix.

cells <- colnames(master_matrix)

sathy_list <- rep("sathyamurthy", each = length(colnames(sathy_full_mat)))

#haring_list <- rep("haring", length(colnames(haring_full_mat))

blum_list <- rep("blum", each = length(colnames(blum_full_mat)))

alkas_list <- rep("alkaslasi", each = length(colnames(alkas_full_mat)))

#study <- c(sathy_list, haring_list, blum_list, alkas_list)
study <- c(sathy_list, blum_list, alkas_list)

metadata <- data.frame(cells, study)
write.csv(metadata, file = "../raw_data_full/metadata.csv")

#seperate by sex from cells using via last underscore using regex expression.
#computationally intensive, need HPC server.
#metadata %>% separate(cells, c("cells", "sex"), sep="_(?=[^_]+$)") 
#write.csv(metadata, file = "../raw_data_full/metadata_split.csv")