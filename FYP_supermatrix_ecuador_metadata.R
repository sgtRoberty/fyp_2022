library(ape)
library(seqinr)
library(readr)
library(dplyr)

AllMitoMetadata <- read_csv("AllMitogenomesMetadata.csv")
supermatrix_ecuador <- read.nexus.data("")

supermatrix_ecuador_metadata <- filter(AllMitoMetadata, 
                                       db_id %in% names(supermatrix_ecuador))
length(names(supermatrix_ecuador)[! names(supermatrix_ecuador) %in% AllMitoMetadata$db_id])

write.fasta(supermatrix_ecuador, names = names(supermatrix_ecuador54dew))
