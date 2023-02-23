
library(Seurat)
library(stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

metadata <- read.csv("../raw_data_full/metadata_split.csv")

metadata$sex <- as.factor(metadata$sex)
metadata$study <- as.factor(metadata$study)
#metadata$cell_type <- as.factor(metadata$cell_type)
#metadata$orig.ident <- metadata$study

sathy_meta <- metadata[metadata$study == "sathyamurthy",]
blum_meta <- metadata[metadata$study == "blum",]
alkaslasi_meta <- metadata[metadata$study == "alkaslasi",]

#add sample identification to metadata.

#sathy
sample_sathy <- str_extract(sathy_meta$cells, "[^_]+")
unique(sample_sathy)

sathy_meta$sample <- as.factor(sample_sathy)

#blum
#they homogenized two samples per replicate.

sample_blum <- str_extract(blum_meta$cells, "[^_]*_[^_]*")
unique(sample_blum)

sample_5_6 <- sample_blum[sample_blum == "5_6"]
sample_5_7 <- sample_blum[sample_blum == "5_7"]
sample_4_4 <- sample_blum[sample_blum == "4_4"]

#sample_5_6[] <- "s56"
#sample_5_7[] <- "s57"
#sample_4_4[] <- "s44"


split_randomize <- function(vector, vector_items, split_length, seed_number){
  #make half (or any given number in "split_length") the vector a different sample
  
  vector_split <- rep(vector_items, each = length(vector)/ split_length, length.out = length(vector))
  
  #randomize vector (with constant seed for reproducibility)
  set.seed(seed_number)
  final_vector <- sample(vector_split)
  return(final_vector)
}

sample_5_6 <- split_randomize(sample_5_6, c("blum1","blum2"), 2, 005)
sample_5_7 <- split_randomize(sample_5_7, c("blum3","blum4"), 2, 006)
sample_4_4 <- split_randomize(sample_4_4, c("blum5","blum6"), 2, 007)

blum_meta$sample <- as.factor(c(sample_5_6, sample_5_7, sample_4_4))

#alkaslasi seperated into 6 males 6 females.

sample_alkas <- str_extract(alkaslasi_meta$cells, "[^_]*_[^_]*")
unique(sample_alkas)

sample_f <- sample_alkas[sample_alkas == "lumb_f"] #lumb_f might not be the correct identifier
sample_m <- sample_alkas[sample_alkas == "lumb_m"] #lumb_m might not be the correct identifier

#sample_f[] <- "slumbf"
#sample_m[] <- "slumbm"

sample_f <- split_randomize(sample_f, c("alkas1", "alkas2", "alkas3", "alkas4", "alkas5", "alkas6"), 6, 001)
sample_m <- split_randomize(sample_m, c("alkas7", "alkas8", "alkas9", "alkas10", "alkas11", "alkas12"), 6, 002)

alkaslasi_meta$sample <- as.factor(c(sample_f, sample_m))
#concat into proper order..
unique(metadata$study)

metadata <- rbind(sathy_meta, blum_meta, alkaslasi_meta)

write.csv(metadata, "../raw_data_full/metadata_final.csv")
