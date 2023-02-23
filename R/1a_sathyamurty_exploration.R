library(Matrix)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

matrix <- read.table(file = "../raw_data_full/sathyamurthy/GSE103892_Expression_Count_Matrix.txt", header = TRUE)

rownames(matrix) <- matrix[,1]
matrix <- matrix[,2:18001]

cells <- colnames(matrix)
genes <- rownames(matrix)

data_matrix <- data.matrix(matrix)

counts_sathyamurthy <- as(data_matrix, "sparseMatrix")
colnames(counts_sathyamurthy) <- cells
rownames(counts_sathyamurthy) <- genes

#Distinguish female and male with this given data from the authors.
#Rota1	F
#Rota2	M
#Rota3	M
#Rota4	F
#Rota5	M

#Form1,6	M
#Form2,7	M
#Form3,8	F
#Form4,9	M
#Form5,10	F

#starts by "M" means male, starts by "F" means female

string <- data.frame(strsplit(cells, split = "_"))
identification <- string[1,]

unique(as.list(identification))

identification[3500]

for (i in 1:length(cells)) {
  if (identification[i] %in% c("m1", "M4", "M5", "Macx", "form1", "form2", "form4", "form6", "form7", "form9", "rotarod2", "rotarod3", "rotarod5")){
    identification[i] <- paste0(cells[i], "_M")
  } else {
    identification[i] <- paste0(cells[i], "_F")
  }
}

colnames(counts_sathyamurthy) <- identification

test <- colnames(counts_sathyamurthy)

write.csv(colnames(counts_sathyamurthy), file = "../raw_data_full/sathyamurthy/cells.csv")
write.csv(rownames(counts_sathyamurthy), file = "../raw_data_full/sathyamurthy/genes.csv")

writeMM(counts_sathyamurthy, file = "../raw_data_full/sathyamurthy/matrix_sathyamurthy.txt")
