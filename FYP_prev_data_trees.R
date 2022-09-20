library(ape)
library(seqinr)
library(readr)
library(dplyr)
library(phytools)

##### Metadata #####
AllMitoMetadata <- read_csv("../data/SITE-100_Mitogenome_Metadata_2022-05-17.csv",
                            col_types = cols(subgenus = "c", 
                                             notes = "c", 
                                             authority = "c"))

unique(substr(AllMitoMetadata$db_id, 1, 4))
# "BGLP" "BIOD" "CCCP" "CDBP" "CERN" "GBDL" "HNSP" "MIZA" "QINL" "SPSO" "SRAA"

BIODMetadata <- filter(AllMitoMetadata,
                       substr(db_id, 1, 4) == "BIOD") # 2995 taxa
table(BIODMetadata$country)
# French Guiana      Malaysia        Panama 
# 427                1311            1104 

metadata_problems <- problems(read_csv("../data/SITE-100_Mitogenome_Metadata_2022-05-17.csv",
                                       col_types = cols(subgenus = "c", 
                                                        notes = "c", 
                                                        authority = "c")))
metadata_problems$file <- NULL
metadata_problems$Sequence_ID <- AllMitoMetadata$db_id[metadata_problems$row - 1]
metadata_problems$Error_Col <- names(AllMitoMetadata[metadata_problems$col])
write.csv(metadata_problems, "../results/Metadata_problems.csv")

##### 5709-taxon mitogenome nucleotide data and tree #####
mito_5709_data <- read.fasta("../data/5709taxa_supermatrix.fasta", 
                             seqtype = "DNA",
                             forceDNAtolower = F) # 5709 taxa, seqlen: 18603
mito_5709_besttree <- read.newick("../data/RAxML_bestTree.result.5709taxa.renamed.tre")

mito_5709_metadata <- dplyr::filter(AllMitoMetadata, 
                                    db_id %in% names(mito_5709_data)) # 5642 taxa
length(names(mito_5709_data)[! names(mito_5709_data) %in% AllMitoMetadata$db_id])
# 67 taxa in the mito_5709_data don't have a match in AllMitoMetadata

table(mito_5709_metadata$family)
table(mito_5709_metadata$country)

##### Malaise trap data and tree #####
global_malaise_data <- read.fasta("../data/GlobalMalaise_aligned.fasta", 
                                  seqtype = "DNA",
                                  forceDNAtolower = F) # 12924 taxa, seqlen: 703
malaise_tree <- read.newick("../data/RAxML_result_renamed.malais_out_all2.tre")

Ntip(malaise_tree) - length(global_malaise_data) # 18633 - 12924 = 5709

MalaiseTrapMetadata <- read_csv("../data/MalaiseTrapMetadata.csv")
names(MalaiseTrapMetadata)[c(1, 2, 3, 4, 5)] <- c("db_id", "taxonomy_notes", 
                                                  "family", "country", "subregion")

table(MalaiseTrapMetadata$family)
table(MalaiseTrapMetadata$country)
length(MalaiseTrapMetadata$db_id) # 12924 taxa

# De-name data and tree
for (n in 1 : length(global_malaise_data)){
  names(global_malaise_data)[n] <- unlist(strsplit(names(global_malaise_data)[n], "\\|"))[1]
}
write.fasta(global_malaise_data, names = names(global_malaise_data),
            file.out = "../data/GlobalMalaise_aligned_noname.fasta")

for (n in 1 : Ntip(malaise_tree)){
  malaise_tree$tip.label[n] <- unlist(strsplit(malaise_tree$tip.label[n], "\\~"))[1]
}
write.tree(malaise_tree, file = "../data/malaise_noname.tre")

##### 5709-taxon mitogenome + Malaise trap data and tree #####
mito_5709_global_malaise_data <- read.fasta("../data/5709taxa_GlobalMalaise_supermatrix.fasta",
                                            seqtype = "DNA",
                                            forceDNAtolower = F) # 18633 taxa, seqlen: 18603

# Branch length estimation by RAxML-ng
brlen_5709taxa_malaise_tree <- read.newick(
  "../data/brlen_5709taxa_malaise.raxml.bestTree") # 18633 taxa
brlen_5709taxa_malaise_lbremoved_tree <- read.newick(
  "../data/brlen_5709taxa_malaise_lbremoved.raxml.bestTree") # 18622 taxa

# Remove 11 exceptionally long branches in Dendroscope
Ntip(brlen_5709taxa_malaise_tree) - Ntip(brlen_5709taxa_malaise_lbremoved_tree) # 11 taxa
removed_longbr_names <- brlen_5709taxa_malaise_tree$tip.label[!brlen_5709taxa_malaise_tree$tip.label %in% brlen_5709taxa_malaise_lbremoved_tree$tip.label]
removed_longbr_names

# 5709-taxon and global Malaise trap metadata
mito_5709_global_malaise_metadata <- bind_rows(mito_5709_metadata, 
                                               MalaiseTrapMetadata) # 18566 taxa = 12924 taxa + 5642 taxa
write_csv(mito_5709_global_malaise_metadata, "../data/5709taxa_GlobalMalaise_Metadata.csv")

paste(brlen_5709taxa_malaise_tree$tip.label[! brlen_5709taxa_malaise_tree$tip.label %in% mito_5709_global_malaise_metadata$db_id], 
      collapse = ',')

##### BIBC data and tree #####
BIBC_6282_data <- read.fasta("../data/BIBC_sequences_2021-04-17.fasta", 
                             seqtype = "DNA",
                             forceDNAtolower = F) # 6282 taxa, seqlen: 418
BIBC_besttree <- read.tree("../data/RAxML_bestTree.bibc_seq_out.tre")

Ntip(bibc_besttree) - length(BIBC_6282_data) # 11991 - 6282 = 5709

BIBC_reps_data <- read.fasta("../data/BIBC_representatives.fasta",
                             seqtype = "DNA",
                             forceDNAtolower = F) # 3580 taxa, seqlen: 418

BIBCMetadata <- read_csv("../data/BIBC_metadata_2020-11-01_TJC_2021-04-02_apv.csv",
                         col_types = cols(subtribe = "c",
                                          specimen_id = "c")) # 7348 taxa
names(BIBCMetadata)[4] <- "db_id"
length(BIBCMetadata$db_id) # 7348 taxa

table(BIBCMetadata$country)
# French Guiana      Malaysia        Panama 
# 2093               3008            2241 

BIBC_6282_metadata <- filter(BIBCMetadata, 
                             db_id %in% names(BIBC_6282_data)) # 6282 taxa
length(BIBCMetadata$db_id[! BIBCMetadata$db_id %in% names(BIBC_6282_data)])
# 1066 taxa in the BIBCMetadata don't have a match in BIBC_6282_data

table(BIBC_6282_metadata$family)
table(BIBC_6282_metadata$country)
# French Guiana      Malaysia        Panama 
# 1480               2805            1997 

BIBC_reps_metadata <- filter(BIBCMetadata, 
                             db_id %in% names(BIBC_reps_data)) # 3580 taxa
length(BIBCMetadata$db_id[! BIBCMetadata$db_id %in% names(BIBC_reps_data)])
# 3768 taxa in the BIBCMetadata don't have a match in BIBC_reps_data
table(BIBC_reps_metadata$family)
table(BIBC_reps_metadata$country)
# French Guiana      Malaysia        Panama 
# 591                1658            1331

##### 5709-taxon mitogenome + BIBC data #####
# 5709-taxon and BIBC metadata
mito_5709_bibc_metadata <- bind_rows(mito_5709_metadata, 
                                     BIBCMetadata) # 12990 taxa = 5642 taxa + 7348 taxa
write_csv(mito_5709_bibc_metadata, "../data/5709taxa_BIBC_Metadata.csv")

##### Ecuador barcodes and OTUs and metadata #####
Ecuador_barcodes <- read.fasta("../data/Ecuador_barcodes_clean.fasta",
                               seqtype = "DNA",
                               forceDNAtolower = F) # 3073 taxa (7 duplicates removed), seqlen: 418
Ecuador_otus <- read.fasta("../data/Ecuador_otus.fasta", 
                           seqtype = "DNA",
                           forceDNAtolower = F) # 2333 taxa, seqlen: 418

Ecuador_barcodes_metadata <- read_csv("../data/Ecuador_database_metadata.csv") # 3671 taxa
colnames(Ecuador_barcodes_metadata) <- chartr(" ", "_", 
                                              names(Ecuador_barcodes_metadata))
colnames(Ecuador_barcodes_metadata) <- tolower(colnames(Ecuador_barcodes_metadata))

Ecuador_barcodes_metadata$db_id <- substr(Ecuador_barcodes_metadata$NHMUK_specimen_barcode,
                                             7, 14)

