library(ape)
library(seqinr)
library(phytools)
library(ggtree)
library(ggplot2)
library(tidytree)
library(readr)
library(dplyr)
library(jpeg)
library(png)
library(cowplot)

##### Metadata #####
AllMitoMetadata <- readr::read_csv("../data/SITE-100_Mitogenome_Metadata_2022-05-17.csv",
                                   col_types = cols(subgenus = "c", 
                                                    notes = "c", 
                                                    authority = "c"))
aa_supermatrix_raw <- read.fasta("../data/5_aa_supermatrix.fasta", 
                                 seqtype = "AA") # 5682 taxa, seqlen: 6372
aa_supermatrix_metadata <- filter(AllMitoMetadata, 
                                  db_id %in% names(aa_supermatrix_raw)) # 5609 taxa
sort(table(aa_supermatrix_metadata$country))
unique(aa_supermatrix_metadata$country)
Malaysia_metadata <- filter(aa_supermatrix_metadata, country == "Malaysia")
Panama_metadata <- filter(aa_supermatrix_metadata, country == "Panama")

##### RAxML-ng output tree #####
raxmlng_intermd_tree <- 
  read.tree("../data/allmito_nostarttree_constr_part.raxml.lastTree.TMP")
raxmlng_intermd_tree_plot <- ggtree(raxmlng_intermd_tree, size = 0.2) %<+% 
  AllMitoMetadata + geom_tree(aes(color=superfamily), size = 0.2) + 
  geom_tippoint(aes(color=superfamily), alpha = 0.6, size = 0.1) + 
  geom_tiplab(size = 0.2, nudge_x = -0.04) +
  theme(legend.position="right")
ggsave("../results/raxmlng_intermd_superfamilytree.pdf", 
       plot = raxmlng_intermd_tree_plot, 
       width = 20, height = 40)

raxmlng_intermd_tree_plot <- ggtree(raxmlng_intermd_tree, size = 0.2) %<+% 
  AllMitoMetadata + geom_tree(aes(color=family), size = 0.2) + 
  geom_tippoint(aes(color=family), alpha = 0.6, size = 0.1) + 
  geom_tiplab(size = 0.2, nudge_x = -0.04) +
  theme(legend.position="right")
ggsave("../results/raxmlng_intermd_familytree.pdf", 
       plot = raxmlng_intermd_tree_plot, 
       width = 20, height = 40)

raxmlng_intermd_tree_RN <- 
  read.newick("../data/allmito_nostarttree_constr_part.raxml.lastTree.TMP_RN.tre")
raxmlng_intermd_tree_RN_plot <- ggtree(raxmlng_intermd_tree_RN, size = 0.1) +
  geom_tiplab(size = 0.2, nudge_x = -0.04)
ggsave("../results/raxmlng_intermd_RN_tree.pdf", 
       plot = raxmlng_intermd_tree_RN_plot, 
       width = 20, height = 50, limitsize = FALSE)

MP_metadata <- rbind(Malaysia_metadata, Panama_metadata)
ultra_allmito_extd <- 
  read.tree("../results/randsample_5682tree/ultra_allmito_5682taxa_tree.tre")
region_coloured_tree_plot <- ggtree(ultra_allmito_extd, 
                                    size = 0.2
                                    ) %<+% 
  MP_metadata + geom_tree(aes(color=country), size = 0.2) + 
  geom_tippoint(aes(color=country), alpha = 0.6, size = 0.1) + 
  geom_tiplab(size = 0.1) +
  scale_colour_manual(values=c('forestgreen','orange'),
                      labels = c("Malaysia", "Panama", "Others"),
                      name = "Country") +
  geom_cladelab(node=1005, label="test label", align=TRUE,  
                  offset = .2, textcolor='red', barcolor='red') +
  geom_cladelab(node=3034, label="another clade", align=TRUE, 
                offset = .2, textcolor='blue', barcolor='blue')+
  theme(legend.position="right")
ggsave("../results/FYP_figs/region_coloured_tree.pdf", 
       plot = region_coloured_tree_plot, 
       width = 20, height = 40)

"../results/mccr_test/ultra_aa_mitogenome_Malaysia_tree.tre"
"../results/mccr_test/ultra_aa_mitogenome_Panama_tree.tre"

##### Superfamily tree, single colour #####

nodelabel_df <- read_csv("../results/FYP_figs/superfamily_nodelabel_df_clean.csv")
nodelabel_df_for_points <- nodelabel_df
nodelabel_df <- filter(nodelabel_df, 
                       ! node_name %in% c("Cucujoidea",
                                          "Byrrhoidea",
                                          "Scirtoidea"))

collapsed_tree_plot <- ggtree(collapsed_superfamily_tree,
                              size = 0.8) +
  geom_point2(aes(subset=node%in%nodelabel_df_for_points$node_number), 
              shape=23, size=1, fill='yellow') + 
  hexpand(0.1) + theme_tree2() + 
  xlab(label = "Time from the oldest root (Millions of years)") +
  theme(axis.title.x = element_text(size = 14, vjust = -1.5, face="bold"), 
        axis.line = element_line(color='black', size = 1),
        axis.text.x = element_text(size=12, vjust = -1))

collapsed_tree_plot <-
  scaleClade(collapsed_tree_plot, 6150, 0.2) %>%
  scaleClade(8911, 0.4) %>%
  scaleClade(7179, 0.2) %>%
  scaleClade(9192, 0.2) %>%
  scaleClade(10598, 0.2) %>%
  scaleClade(8331, 0.2) %>%
  scaleClade(11057, 20) %>% # Rhinorhipoidea
  scaleClade(10594, 15) %>%
  scaleClade(11100, 15) %>%
  scaleClade(11093, 15) %>%
  scaleClade(10592, 15) %>% # Derodontoidea
  scaleClade(5666, 25) %>% # Raphidioptera
  scaleClade(11325, 15) %>% # Neuroptera
  scaleClade(11324, 15) %>% # Strepsiptera
  scaleClade(8326, 9) + # Lymexyloidea
  geom_cladelab(data = nodelabel_df, 
                mapping = aes(node = node_number, 
                              label = node_name),
                fontcolor = "steelblue",barcolor = "steelblue",
                fontsize = 4.5, barsize = 3, align = T) +
  geom_cladelab(node = 5999, label = "Cucujoidea",
                fontsize = 4.5, barsize = 3, vjust = 2,
                textcolor="orange", barcolor="orange") +
  geom_cladelab(node = 8162, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor="orange") +
  geom_cladelab(node = 8222, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor="orange") +
  geom_cladelab(node = 8321, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor="orange") +
  geom_cladelab(node = 10928, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor="orange") +
  geom_cladelab(node = 10931, label = "Byrrhoidea", 
                fontsize = 4.5, barsize = 3, vjust = 0,
                textcolor="orange", barcolor="orange") +
  geom_cladelab(node = 11058, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor = "orange") +
  geom_cladelab(node = 11067, label = "Scirtoidea", 
                fontsize = 4.5, barsize = 3, vjust = 0.2,
                textcolor="orange", barcolor="orange") +
  theme(legend.position = "none", 
        plot.margin = grid::unit(c(0,10,6,10), "mm"))


collapsed_tree_plot <- collapsed_tree_plot %>%
  collapse(node= 5666 , 'mixed', fill = "steelblue") %>%
  collapse(node= 5676 , 'mixed', fill = "steelblue") %>%
  collapse(node= 5755 , 'mixed', fill = "steelblue") %>%
  collapse(node= 5999 , 'mixed', fill = "steelblue") %>% # Cucujoidea
  collapse(node= 6150 , 'mixed', fill = "steelblue") %>%
  collapse(node= 7179 , 'mixed', fill = "steelblue") %>%
  collapse(node= 8162 , 'mixed', fill = "steelblue") %>% # Cucujoidea
  collapse(node= 8222 , 'mixed', fill = "steelblue") %>% # Cucujoidea
  collapse(node= 8321 , 'mixed', fill = "steelblue") %>% # Cucujoidea
  collapse(node= 8326 , 'mixed', fill = "steelblue") %>%
  collapse(node= 8331 , 'mixed', fill = "steelblue") %>%
  collapse(node= 8784 , 'mixed', fill = "steelblue") %>%
  collapse(node= 8911 , 'mixed', fill = "steelblue") %>%
  collapse(node= 9192 , 'mixed', fill = "steelblue") %>%
  collapse(node= 10408 , 'mixed', fill = "steelblue") %>%
  collapse(node= 10506 , 'mixed', fill = "steelblue") %>%
  collapse(node= 10592 , 'mixed', fill = "steelblue") %>%
  collapse(node= 10594 , 'mixed', fill = "steelblue") %>%
  collapse(node= 10598 , 'mixed', fill = "steelblue") %>%
  collapse(node= 10928 , 'mixed', fill = 'steelblue') %>%
  collapse(node= 10931 , 'mixed', fill = 'steelblue') %>%
  collapse(node= 11021 , 'mixed', fill = "steelblue") %>%
  collapse(node= 11057 , 'mixed', fill = "steelblue") %>%
  collapse(node= 11058 , 'mixed', fill = 'steelblue') %>%
  collapse(node= 11067 , 'mixed', fill = 'steelblue') %>%
  collapse(node= 11093 , 'mixed', fill = "steelblue") %>%
  collapse(node= 11100 , 'mixed', fill = "steelblue") %>%
  collapse(node= 11103 , 'mixed', fill = "steelblue") %>%
  collapse(node= 11149 , 'mixed', fill = "steelblue") %>%
  collapse(node= 11324 , 'mixed', fill = "steelblue") %>%
  collapse(node= 11325 , 'mixed', fill = "steelblue") +
  scale_fill_discrete(aes(colour = node_name))

ggsave("../results/FYP_figs/test_collapsed_blue.pdf", 
       plot = collapsed_tree_plot, 
       width = 20, height = 20)



##### family, collapsed tree #####
collapsed_tree <- 
  read.tree("../results/randsample_5682tree/ultra_allmito_5682taxa_tree_RN_collapsed.tree")
sort(unique(collapsed_tree$node.label))

# Node labels and node numbers
nodelabel_list <- collapsed_tree$node.label
node_number_list <- c()
node_name_list <- c()
for (i in 1:length(nodelabel_list)){
  if (nodelabel_list[i] != ""){
    node_number_list <- append(node_number_list, i)
    node_name_list <- append(node_name_list, nodelabel_list[i])
  }
}
nodelabel_df <- data.frame(c(node_number_list), c(node_name_list))
colnames(nodelabel_df) <- c("node_number", "node_name")
nodelabel_df$node_number <- nodelabel_df$node_number + 5622

# write.csv(nodelabel_df, "../results/FYP_figs/nodelabel_df.csv")

for (number in nodelabel_df$node_number){
  print(paste("collapse(node=", number, ", 'mixed') %>%"), quote=F)
}

nodelabel_df <- read_csv("../results/FYP_figs/nodelabel_df.csv")

collapsed_tree_plot <- ggtree(collapsed_tree, 
                              size = 0.8) +
  geom_point2(aes(subset=node%in%nodelabel_df$node_number), 
              shape=23, size=1, fill='darkgreen')+
  geom_tiplab(size = 0.2, nudge_x = -0.5)

collapsed_tree_plot <-
  scaleClade(collapsed_tree_plot, 6105, 0.2) %>%
  scaleClade(9223, 0.2) %>%
  scaleClade(7229, 0.2) + 
  geom_cladelab(data = nodelabel_df, 
                mapping = aes(node = node_number, 
                              label = node_name,
                              colour = node_name), 
                fontsize = 3) +
  theme(legend.position = "none")

collapsed_tree_plot <- collapsed_tree_plot %>%
  collapse(node= 5624 , 'mixed') %>%
  collapse(node= 5635 , 'mixed') %>%
  collapse(node= 5658 , 'mixed') %>%
  collapse(node= 5698 , 'mixed') %>%
  collapse(node= 5703 , 'mixed') %>%
  collapse(node= 5715 , 'mixed') %>%
  collapse(node= 5722 , 'mixed') %>%
  collapse(node= 5737 , 'mixed') %>%
  collapse(node= 5742 , 'mixed') %>%
  collapse(node= 5759 , 'mixed') %>%
  collapse(node= 5780 , 'mixed') %>%
  collapse(node= 5896 , 'mixed') %>%
  collapse(node= 5943 , 'mixed') %>%
  collapse(node= 5956 , 'mixed') %>%
  collapse(node= 5957 , 'mixed') %>%
  collapse(node= 6089 , 'mixed') %>%
  collapse(node= 6105 , 'mixed') %>%
  collapse(node= 6900 , 'mixed') %>%
  collapse(node= 6908 , 'mixed') %>%
  collapse(node= 6913 , 'mixed') %>%
  collapse(node= 6919 , 'mixed') %>%
  collapse(node= 6922 , 'mixed') %>%
  collapse(node= 7029 , 'mixed') %>%
  collapse(node= 7096 , 'mixed') %>%
  collapse(node= 7100 , 'mixed') %>%
  collapse(node= 7137 , 'mixed') %>%
  collapse(node= 7141 , 'mixed') %>%
  collapse(node= 7177 , 'mixed') %>%
  collapse(node= 7229 , 'mixed') %>%
  collapse(node= 8018 , 'mixed') %>%
  collapse(node= 8036 , 'mixed') %>%
  collapse(node= 8049 , 'mixed') %>%
  collapse(node= 8074 , 'mixed') %>%
  collapse(node= 8112 , 'mixed') %>%
  collapse(node= 8175 , 'mixed') %>%
  collapse(node= 8220 , 'mixed') %>%
  collapse(node= 8240 , 'mixed') %>%
  collapse(node= 8264 , 'mixed') %>%
  collapse(node= 8265 , 'mixed') %>%
  collapse(node= 8269 , 'mixed') %>%
  collapse(node= 8274 , 'mixed') %>%
  collapse(node= 8283 , 'mixed') %>%
  collapse(node= 8441 , 'mixed') %>%
  collapse(node= 8445 , 'mixed') %>%
  collapse(node= 8457 , 'mixed') %>%
  collapse(node= 8469 , 'mixed') %>%
  collapse(node= 8484 , 'mixed') %>%
  collapse(node= 8488 , 'mixed') %>%
  collapse(node= 8502 , 'mixed') %>%
  collapse(node= 8510 , 'mixed') %>%
  collapse(node= 8513 , 'mixed') %>%
  collapse(node= 8526 , 'mixed') %>%
  collapse(node= 8529 , 'mixed') %>%
  collapse(node= 8538 , 'mixed') %>%
  collapse(node= 8555 , 'mixed') %>%
  collapse(node= 8560 , 'mixed') %>%
  collapse(node= 8561 , 'mixed') %>%
  collapse(node= 8606 , 'mixed') %>%
  collapse(node= 8607 , 'mixed') %>%
  collapse(node= 8713 , 'mixed') %>%
  collapse(node= 8747 , 'mixed') %>%
  collapse(node= 8792 , 'mixed') %>%
  collapse(node= 8800 , 'mixed') %>%
  collapse(node= 8835 , 'mixed') %>%
  collapse(node= 8838 , 'mixed') %>%
  collapse(node= 8870 , 'mixed') %>%
  collapse(node= 8873 , 'mixed') %>%
  collapse(node= 8875 , 'mixed') %>%
  collapse(node= 8890 , 'mixed') %>%
  collapse(node= 9110 , 'mixed') %>%
  collapse(node= 9113 , 'mixed') %>%
  collapse(node= 9193 , 'mixed') %>%
  collapse(node= 9194 , 'mixed') %>%
  collapse(node= 9223 , 'mixed') %>%
  collapse(node= 10325 , 'mixed') %>%
  collapse(node= 10425 , 'mixed') %>%
  collapse(node= 10502 , 'mixed') %>%
  collapse(node= 10505 , 'mixed') %>%
  collapse(node= 10506 , 'mixed') %>%
  collapse(node= 10509 , 'mixed') %>%
  collapse(node= 10513 , 'mixed') %>%
  collapse(node= 10518 , 'mixed') %>%
  collapse(node= 10524 , 'mixed') %>%
  collapse(node= 10530 , 'mixed') %>%
  collapse(node= 10574 , 'mixed') %>%
  collapse(node= 10644 , 'mixed') %>%
  collapse(node= 10748 , 'mixed') %>%
  collapse(node= 10755 , 'mixed') %>%
  collapse(node= 10761 , 'mixed') %>%
  collapse(node= 10808 , 'mixed') %>%
  collapse(node= 10842 , 'mixed') %>%
  collapse(node= 10848 , 'mixed') %>%
  collapse(node= 10850 , 'mixed') %>%
  collapse(node= 10851 , 'mixed') %>%
  collapse(node= 10854 , 'mixed') %>%
  collapse(node= 10871 , 'mixed') %>%
  collapse(node= 10878 , 'mixed') %>%
  collapse(node= 10921 , 'mixed') %>%
  collapse(node= 10926 , 'mixed') %>%
  collapse(node= 10934 , 'mixed') %>%
  collapse(node= 10972 , 'mixed') %>%
  collapse(node= 10973 , 'mixed') %>%
  collapse(node= 10979 , 'mixed') %>%
  collapse(node= 11007 , 'mixed') %>%
  collapse(node= 11014 , 'mixed') %>%
  collapse(node= 11019 , 'mixed') %>%
  collapse(node= 11025 , 'mixed') %>%
  collapse(node= 11058 , 'mixed') %>%
  collapse(node= 11059 , 'mixed') %>%
  collapse(node= 11061 , 'mixed') %>%
  collapse(node= 11076 , 'mixed') %>%
  collapse(node= 11235 , 'mixed') %>%
  collapse(node= 11239 , 'mixed') %>%
  collapse(node= 11241 , 'mixed') 

ggsave("../results/FYP_figs/test_collapsed.pdf", 
       plot = collapsed_tree_plot, 
       width = 20, height = 40)

##### superfamily, collapsed tree #####
collapsed_superfamily_tree <- 
  read.tree("../results/randsample_5682tree/ultra_allmito_5682taxa_tree_RN_superfamily_collapsed.tree")

# Node labels and node numbers
nodelabel_list <- collapsed_superfamily_tree$node.label
node_number_list <- c()
node_name_list <- c()
for (i in 1:length(nodelabel_list)){
  if (nodelabel_list[i] != ""){
    node_number_list <- append(node_number_list, i)
    node_name_list <- append(node_name_list, nodelabel_list[i])
  }
}
nodelabel_df <- data.frame(c(node_number_list), c(node_name_list))
colnames(nodelabel_df) <- c("node_number", "node_name")
nodelabel_df$node_number <- nodelabel_df$node_number + 5664

write.csv(nodelabel_df, "../results/FYP_figs/superfamily_nodelabel_df.csv")

for (number in nodelabel_df$node_number){
  print(paste("collapse(node=", number, ", 'mixed') %>%"), quote=F)
}


nodelabel_df <- read_csv("../results/FYP_figs/superfamily_nodelabel_df_clean.csv")
nodelabel_df_for_points <- nodelabel_df
nodelabel_df <- filter(nodelabel_df, 
                       ! node_name %in% c("Cucujoidea",
                                          "Byrrhoidea",
                                          "Scirtoidea"))

collapsed_tree_plot <- ggtree(collapsed_superfamily_tree,
                              size = 0.8) +
  geom_point2(aes(subset=node%in%nodelabel_df_for_points$node_number), 
              shape=23, size=1, fill='yellow') + 
  hexpand(0.2) + theme_tree2() + 
  xlab(label = "Time from the oldest root (Millions of years)") +
  theme(axis.title.x = element_text(size = 18, vjust = -1.5, face="bold"), 
        axis.line = element_line(color='black', size = 1),
        axis.text.x = element_text(size=16, vjust = -1))

collapsed_tree_plot <-
  scaleClade(collapsed_tree_plot, 6150, 0.3) %>% # Chrysomeloidea
  scaleClade(8911, 0.5) %>%
  scaleClade(7179, 0.3) %>%
  scaleClade(9192, 0.3) %>%
  scaleClade(10598, 0.3) %>%
  scaleClade(8331, 0.3) %>%
  scaleClade(11057, 20) %>% # Rhinorhipoidea
  scaleClade(10594, 15) %>%
  scaleClade(11100, 15) %>%
  scaleClade(11093, 15) %>%
  scaleClade(10592, 15) %>% # Derodontoidea
  scaleClade(5666, 25) %>% # Raphidioptera
  scaleClade(11325, 15) %>% # Neuroptera
  scaleClade(11324, 15) %>% # Strepsiptera
  scaleClade(8326, 10) # Lymexyloidea 

collapsed_tree_plot <-
  collapsed_tree_plot + 
  geom_cladelab(data = nodelabel_df, 
                mapping = aes(node = node_number, 
                              label = node_name,
                              colour = node_name),
                fontsize = 4.5, barsize = 3, align = T, fontface = 2) +
  geom_cladelab(node = 5999, label = "Cucujoidea",
                fontsize = 4.5, barsize = 3, vjust = 2, extend = 5,
                textcolor='#00b7e8', barcolor='#00b7e8', fontface = 2) +
  geom_cladelab(node = 8162, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor='#00b7e8') +
  geom_cladelab(node = 8222, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor='#00b7e8', extend = 1.1) +
  geom_cladelab(node = 8321, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor='#00b7e8') +
  geom_cladelab(node = 10928, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor='#ff6c92', extend = 0.6) +
  geom_cladelab(node = 10931, label = "Byrrhoidea", 
                fontsize = 4.5, barsize = 3, vjust = 0, extend = 0.5,
                textcolor='#ff6c92', barcolor='#ff6c92', fontface = 2) +
  geom_cladelab(node = 11058, label = "", 
                fontsize = 4.5, barsize = 3,
                barcolor = '#00aefa', extend = 0.6) +
  geom_cladelab(node = 11067, label = "Scirtoidea", 
                fontsize = 4.5, barsize = 3, vjust = 0.2, extend = 0.5,
                textcolor="#00aefa", barcolor="#00aefa", fontface = 2) +
  theme(legend.position = "none", 
        plot.margin = grid::unit(c(0,10,6,10), "mm"))

collapsed_tree_plot <- collapsed_tree_plot %>%
  collapse(node= 5666 , 'mixed', fill = "#be80ff") %>%
  collapse(node= 5676 , 'mixed', fill = "#f8766d") %>%
  collapse(node= 5755 , 'mixed', fill = "#97a900") %>%
  collapse(node= 5999 , 'mixed', fill = "#00b7e8") %>% # Cucujoidea
  collapse(node= 6150 , 'mixed', fill = "#ca9700") %>%
  collapse(node= 7179 , 'mixed', fill = "#71b000") %>%
  collapse(node= 8162 , 'mixed', fill = "#00b7e8") %>% # Cucujoidea
  collapse(node= 8222 , 'mixed', fill = "#00b7e8") %>% # Cucujoidea
  collapse(node= 8321 , 'mixed', fill = "#00b7e8") %>% # Cucujoidea
  collapse(node= 8326 , 'mixed', fill = "#3da1ff") %>%
  collapse(node= 8331 , 'mixed', fill = "#ff6c92") %>%
  collapse(node= 8784 , 'mixed', fill = "#b3a000") %>%
  collapse(node= 8911 , 'mixed', fill = "#f265e7") %>%
  collapse(node= 9192 , 'mixed', fill = "#fe61cf") %>%
  collapse(node= 10408 , 'mixed', fill = "#00b7e8") %>%
  collapse(node= 10506 , 'mixed', fill = "#00aefa") %>%
  collapse(node= 10592 , 'mixed', fill = "#00bb4b") %>%
  collapse(node= 10594 , 'mixed', fill = "#2fb600") %>%
  collapse(node= 10598 , 'mixed', fill = "#00c098") %>%
  collapse(node= 10928 , 'mixed', fill = '#ff6c92') %>%
  collapse(node= 10931 , 'mixed', fill = '#ff6c92') %>%
  collapse(node= 11021 , 'mixed', fill = "#ec823c") %>%
  collapse(node= 11057 , 'mixed', fill = "#de71f9") %>%
  collapse(node= 11058 , 'mixed', fill = '#00aefa') %>%
  collapse(node= 11067 , 'mixed', fill = '#00aefa') %>%
  collapse(node= 11093 , 'mixed', fill = "#00c0b7") %>%
  collapse(node= 11100 , 'mixed', fill = "#00bdd1") %>%
  collapse(node= 11103 , 'mixed', fill = "#00bf76") %>%
  collapse(node= 11149 , 'mixed', fill = "#dd8d00") %>%
  collapse(node= 11324 , 'mixed', fill = "#ff64b3") %>%
  collapse(node= 11325 , 'mixed', fill = "#8f91ff") +
  scale_fill_discrete(aes(colour = node_name)) 

Chrysomeloidea <- readPNG("../results/FYP_figs/Chrysomeloidea_4.png")
Chrysomeloidea <- rphylopic:::recolor_phylopic(Chrysomeloidea, 
                                               color = "#ca9700", 
                                               alpha = 1)
Chrysomeloidea_2 <- readPNG("../results/FYP_figs/Chrysomeloidea_2.png")
Chrysomeloidea_2 <- rphylopic:::recolor_phylopic(Chrysomeloidea_2, 
                                               color = "#ca9700", 
                                               alpha = 1)
Curculionoidea_1 <- readPNG("../results/FYP_figs/Curculionoidea_1.png")
Curculionoidea_1 <- rphylopic:::recolor_phylopic(Curculionoidea_1, 
                                               color = "#71b000", 
                                               alpha = 1)
Curculionoidea_2 <- readPNG("../results/FYP_figs/Curculionoidea_2.png")
Curculionoidea_2 <- rphylopic:::recolor_phylopic(Curculionoidea_2, 
                                                 color = "#71b000", 
                                                 alpha = 1)
Cucujoidea <- readPNG("../results/FYP_figs/Cucujoidea.png")
Cucujoidea <- rphylopic:::recolor_phylopic(Cucujoidea, 
                                           color = "#00b7e8", 
                                           alpha = 1)
Tenebrionoidea <- readPNG("../results/FYP_figs/Tenebrionoidea.png")
Tenebrionoidea <- rphylopic:::recolor_phylopic(Tenebrionoidea, 
                                           color = "#ff6c92", 
                                           alpha = 1)

Coccinelloidea <- readPNG("../results/FYP_figs/Coccinelloidea.png")
Coccinelloidea <- rphylopic:::recolor_phylopic(Coccinelloidea, 
                                               color = "#97a900", 
                                               alpha = 1)
Staphylinoidea <- readPNG("../results/FYP_figs/Staphylinoidea.png")
Staphylinoidea <- rphylopic:::recolor_phylopic(Staphylinoidea, 
                                               color = "#fe61cf", 
                                               alpha = 1)
Scarabaeoidea <- readPNG("../results/FYP_figs/Scarabaeoidea.png")
Scarabaeoidea <- rphylopic:::recolor_phylopic(Scarabaeoidea, 
                                              color = "#f265e7", 
                                              alpha = 1)
Scarabaeoidea_3 <- readPNG("../results/FYP_figs/Scarabaeoidea_3.png")
Scarabaeoidea_3 <- rphylopic:::recolor_phylopic(Scarabaeoidea_3, 
                                              color = "#f265e7", 
                                              alpha = 1)
Hydrophiloidea <- readPNG("../results/FYP_figs/Hydrophiloidea.png")
Hydrophiloidea <- rphylopic:::recolor_phylopic(Hydrophiloidea, 
                                                color = "#00aefa", 
                                                alpha = 1)
Elateroidea_2 <- readPNG("../results/FYP_figs/Elateroidea_2.png")
Elateroidea_2 <- rphylopic:::recolor_phylopic(Elateroidea_2, 
                                              color = "#00c098", 
                                              alpha = 1)
Byrrhoidea_4 <- readPNG("../results/FYP_figs/Byrrhoidea_4.png")
Byrrhoidea_4 <- rphylopic:::recolor_phylopic(Byrrhoidea_4, 
                                              color = "#ff6c92", 
                                              alpha = 1)

Buprestoidea <- readPNG("../results/FYP_figs/Buprestoidea.png")
Buprestoidea <- rphylopic:::recolor_phylopic(Buprestoidea, 
                                             color = "#ec823c", 
                                             alpha = 1)

Caraboidea <- readPNG("../results/FYP_figs/Caraboidea.png")
Caraboidea <- rphylopic:::recolor_phylopic(Caraboidea, 
                                           color = "#dd8d00", 
                                           alpha = 1)
Dytiscoidea <- readPNG("../results/FYP_figs/Dytiscoidea.png")
Dytiscoidea <- rphylopic:::recolor_phylopic(Dytiscoidea, 
                                            color = "#00bf76", 
                                            alpha = 1)

Gyrinoidea <- readPNG("../results/FYP_figs/Gyrinus_substriatus_Stephens,_1828.png")
Gyrinoidea <- rphylopic:::recolor_phylopic(Gyrinoidea, 
                                            color = "#00c0b7", 
                                            alpha = 1)
Cleroidea <- readPNG("../results/FYP_figs/Tarsostenus_univittatus-Curtis_rendered.png")
Cleroidea <- rphylopic:::recolor_phylopic(Cleroidea, 
                                          color = "#b3a000", 
                                          alpha = 1)
Strepsiptera <- readPNG("../results/FYP_figs/Strepsiptera.png")
Strepsiptera <- rphylopic:::recolor_phylopic(Strepsiptera, 
                                          color = "#ff64b3", 
                                          alpha = 1)

collapsed_tree_plot_image <- 
ggdraw() +
  draw_plot(collapsed_tree_plot) +  
  draw_image(Chrysomeloidea,  
             x = 0.42, y = 0.477, scale = .06) +  
  draw_image(Chrysomeloidea_2,  
             x = 0.425, y = 0.428, scale = .08) +
  draw_image(Curculionoidea_1,  
             x = 0.43, y = 0.37, scale = .08) +
  draw_image(Curculionoidea_2,  
             x = 0.43, y = 0.33, scale = .08) +
  draw_image(Cucujoidea,  
             x = 0.42, y = 0.255, scale = .08) +
  draw_image(Tenebrionoidea,  
             x = 0.435, y = 0.181, scale = .08) +
  draw_image(Coccinelloidea,  
             x = 0.43, y = 0.104, scale = .06) +
  draw_image(Cleroidea,  
             x = 0.43, y = 0.047, scale = .1) +
  draw_image(Staphylinoidea,  
             x = 0.425, y = -0.06, scale = .08) +
  draw_image(Scarabaeoidea,  
             x = 0.438, y = -0.12, scale = .06) +
  draw_image(Scarabaeoidea_3,  
             x = 0.43, y = -0.155, scale = .095) +
  draw_image(Hydrophiloidea,  
             x = 0.43, y = -0.208, scale = .065) +
  draw_image(Elateroidea_2,  
             x = 0.42, y = -0.242, scale = .065) +
  draw_image(Byrrhoidea_4,  
             x = 0.425, y = -0.271, scale = .055) +
  draw_image(Buprestoidea,  
             x = 0.43, y = -0.295, scale = .06) +
  draw_image(Caraboidea,  
             x = 0.43, y = -0.3575, scale = .075) +
  draw_image(Dytiscoidea,  
             x = 0.43, y = -0.393, scale = .056) +
  draw_image(Gyrinoidea,  
             x = 0.43, y = -0.426, scale = .06) +
  draw_image(Strepsiptera,  
             x = 0.43, y = -0.454, scale = .035) +
  draw_line(x = c(0.913, 0.867),
            y = c(0.0465, 0.056),
            color = "#ff64b3", size = 0.5)

ggsave("../results/FYP_figs/coleoptera_superfamily_collapsed_tree.pdf", 
       plot = collapsed_tree_plot_image, 
       width = 20, height = 25)

ggsave("../results/FYP_figs/coleoptera_superfamily_collapsed_tree.png", 
       plot = collapsed_tree_plot_image, 
       width = 20, height = 25)


##### superfamily, collapsed tree BW #####
collapsed_superfamily_tree <- 
  read.tree("../results/randsample_5682tree/ultra_allmito_5682taxa_tree_RN_superfamily_collapsed.tree")

# Node labels and node numbers
nodelabel_list <- collapsed_superfamily_tree$node.label
node_number_list <- c()
node_name_list <- c()
for (i in 1:length(nodelabel_list)){
  if (nodelabel_list[i] != ""){
    node_number_list <- append(node_number_list, i)
    node_name_list <- append(node_name_list, nodelabel_list[i])
  }
}
nodelabel_df <- data.frame(c(node_number_list), c(node_name_list))
colnames(nodelabel_df) <- c("node_number", "node_name")
nodelabel_df$node_number <- nodelabel_df$node_number + 5664

write.csv(nodelabel_df, "../results/FYP_figs/superfamily_nodelabel_df.csv")

for (number in nodelabel_df$node_number){
  print(paste("collapse(node=", number, ", 'mixed') %>%"), quote=F)
}


nodelabel_df <- read_csv("../results/FYP_figs/superfamily_nodelabel_df_clean.csv")
nodelabel_df_for_points <- nodelabel_df
nodelabel_df <- filter(nodelabel_df, 
                       ! node_name %in% c("Cucujoidea",
                                          "Byrrhoidea",
                                          "Scirtoidea"))

collapsed_tree_plot <- ggtree(collapsed_superfamily_tree,
                              size = 0.8) +
  geom_point2(aes(subset=node%in%nodelabel_df_for_points$node_number), 
              shape=23, size=1, fill='yellow') + 
  hexpand(0.2) + theme_tree2() + 
  xlab(label = "Time from the oldest root (Millions of years)") +
  theme(axis.title.x = element_text(size = 20, vjust = -1.5, face="bold"), 
        axis.line = element_line(color='black', size = 1),
        axis.text.x = element_text(size=18, vjust = -1))

collapsed_tree_plot <-
  scaleClade(collapsed_tree_plot, 6150, 0.3) %>% # Chrysomeloidea
  scaleClade(8911, 0.5) %>%
  scaleClade(7179, 0.3) %>%
  scaleClade(9192, 0.3) %>%
  scaleClade(10598, 0.3) %>%
  scaleClade(8331, 0.3) %>%
  scaleClade(11057, 20) %>% # Rhinorhipoidea
  scaleClade(10594, 15) %>%
  scaleClade(11100, 15) %>%
  scaleClade(11093, 15) %>%
  scaleClade(10592, 15) %>% # Derodontoidea
  scaleClade(5666, 25) %>% # Raphidioptera
  scaleClade(11325, 15) %>% # Neuroptera
  scaleClade(11324, 15) %>% # Strepsiptera
  scaleClade(8326, 10) # Lymexyloidea 

collapsed_tree_plot <-
  collapsed_tree_plot + 
  geom_cladelab(data = nodelabel_df, 
                mapping = aes(node = node_number, 
                              label = node_name),
                fontsize = 5.5, barsize = 3, align = T, fontface = 2) +
  geom_cladelab(node = 5999, label = "Cucujoidea",
                fontsize = 5.5, barsize = 3, vjust = 2, extend = 5,
                fontface = 2) +
  geom_cladelab(node = 8162, label = "", 
                fontsize = 5.5, barsize = 3) +
  geom_cladelab(node = 8222, label = "", 
                fontsize = 5.5, barsize = 3,
                extend = 1.1) +
  geom_cladelab(node = 8321, label = "", 
                fontsize = 5.5, barsize = 3)+
  geom_cladelab(node = 10928, label = "", 
                fontsize = 5.5, barsize = 3,
                extend = 0.55) +
  geom_cladelab(node = 10931, label = "Byrrhoidea", 
                fontsize = 5.5, barsize = 3, vjust = 0, extend = 0.5,
                fontface = 2) +
  geom_cladelab(node = 11058, label = "", 
                fontsize = 5.5, barsize = 3,
                extend = 0.5) +
  geom_cladelab(node = 11067, label = "Scirtoidea", 
                fontsize = 5.5, barsize = 3, vjust = 0.2, extend = 0.55,
                fontface = 2) +
  theme(legend.position = "none", 
        plot.margin = grid::unit(c(0,10,6,10), "mm")) + scale_fill_brewer()

collapsed_tree_plot <- collapsed_tree_plot %>%
  collapse(node= 5666 , 'mixed', fill = "#333333") %>% 
  collapse(node= 5676 , 'mixed', fill = "#919191") %>% 
  collapse(node= 5755 , 'mixed', fill = "#919191") %>%
  collapse(node= 5999 , 'mixed', fill = "#333333") %>% # Cucujoidea
  collapse(node= 6150 , 'mixed', fill = "#333333") %>%
  collapse(node= 7179 , 'mixed', fill = "#919191") %>%
  collapse(node= 8162 , 'mixed', fill = "#333333") %>% # Cucujoidea
  collapse(node= 8222 , 'mixed', fill = "#333333") %>% # Cucujoidea
  collapse(node= 8321 , 'mixed', fill = "#333333") %>% # Cucujoidea
  collapse(node= 8326 , 'mixed', fill = "#333333") %>%
  collapse(node= 8331 , 'mixed', fill = "#cccccc") %>%
  collapse(node= 8784 , 'mixed', fill = "#333333") %>%
  collapse(node= 8911 , 'mixed', fill = "#919191") %>%
  collapse(node= 9192 , 'mixed', fill = "#cccccc") %>%
  collapse(node= 10408 , 'mixed', fill = "#333333") %>%
  collapse(node= 10506 , 'mixed', fill = "#cccccc") %>%
  collapse(node= 10592 , 'mixed', fill = "#333333") %>%
  collapse(node= 10594 , 'mixed', fill = "#333333") %>%
  collapse(node= 10598 , 'mixed', fill = "#919191") %>%
  collapse(node= 10928 , 'mixed', fill = '#333333') %>% # Byrrhoidea
  collapse(node= 10931 , 'mixed', fill = '#333333') %>% # Byrrhoidea
  collapse(node= 11021 , 'mixed', fill = "#cccccc") %>%
  collapse(node= 11057 , 'mixed', fill = "#919191") %>%
  collapse(node= 11058 , 'mixed', fill = '#cccccc') %>% # Scirtoidea
  collapse(node= 11067 , 'mixed', fill = '#cccccc') %>% # Scirtoidea
  collapse(node= 11093 , 'mixed', fill = "#333333") %>%
  collapse(node= 11100 , 'mixed', fill = "#919191") %>%
  collapse(node= 11103 , 'mixed', fill = "#cccccc") %>%
  collapse(node= 11149 , 'mixed', fill = "#333333") %>%
  collapse(node= 11324 , 'mixed', fill = "#919191") %>%
  collapse(node= 11325 , 'mixed', fill = "#cccccc")  


collapsed_tree_plot_image <- 
  ggdraw() +
  draw_plot(collapsed_tree_plot) +  
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Chrysomeloidea_4.png"), 
                                          color = "#333333", 
                                          alpha = 1),  
             x = 0.45, y = 0.466, scale = .06) +  
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Chrysomeloidea_2.png"),
                                          color = "#333333", 
                                          alpha = 1), 
             x = 0.453, y = 0.434, scale = .08) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Curculionoidea_1.png"),
                                          color = "#919191", 
                                          alpha = 1),  
             x = 0.45, y = 0.367, scale = .08) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Curculionoidea_2.png"),
                                          color = "#919191", 
                                          alpha = 1),  
             x = 0.45, y = 0.335, scale = .08) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Cucujoidea.png"),
                                          color = "#333333", 
                                          alpha = 1),  
             x = 0.44, y = 0.2535, scale = .08) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Tenebrionoidea.png"),
                                          color = "#cccccc", 
                                          alpha = 1),  
             x = 0.455, y = 0.181, scale = .075) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Coccinelloidea.png"),
                                          color = "#919191", 
                                          alpha = 1),  
             x = 0.453, y = 0.103, scale = .06) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Tarsostenus_univittatus-Curtis_rendered.png"),
                                          color = "#333333",
                                          alpha = 1),  
             x = 0.44, y = 0.045, scale = .13) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Staphylinoidea.png"),
                                          color = "#cccccc",
                                          alpha = 1),  
             x = 0.455, y = -0.063, scale = .08) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Scarabaeoidea.png"),
                                          color = "#919191",
                                          alpha = 1),
             x = 0.46, y = -0.132, scale = .06) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Scarabaeoidea_3.png"),
                                          color = "#919191",
                                          alpha = 1),  
             x = 0.452, y = -0.157, scale = .095) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Hydrophiloidea.png"),
                                          color = "#cccccc",
                                          alpha = 1),  
             x = 0.455, y = -0.21, scale = .065) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Elateroidea_2.png"),
                                          color = "#919191",
                                          alpha = 1),  
             x = 0.43, y = -0.2445, scale = .065) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Byrrhoidea_4.png"),
                                          color = "#333333",
                                          alpha = 1),  
             x = 0.435, y = -0.2735, scale = .055) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Buprestoidea.png"),
                                          color = "#cccccc",   
                                          alpha = 1),  
             x = 0.44, y = -0.298, scale = .06) + 
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Caraboidea.png"),
                                          color = "#333333",   
                                          alpha = 1),  
             x = 0.445, y = -0.36, scale = .075) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Dytiscoidea.png"),
                                          color = "#cccccc",   
                                          alpha = 1),  
             x = 0.44, y = -0.395, scale = .056) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Gyrinus_substriatus_Stephens,_1828.png"),
                                          color = "#333333",   
                                          alpha = 1),  
             x = 0.44, y = -0.43, scale = .06) +
  draw_image(rphylopic:::recolor_phylopic(readPNG("../results/FYP_figs/Strepsiptera.png"),
                                          color = "#919191",   
                                          alpha = 1),  
             x = 0.445, y = -0.4565, scale = .05) +
  draw_line(x = c(0.9195, 0.892),
            y = c(0.044, 0.054),
            size = 0.5, color = "#919191")

ggsave("../results/FYP_figs/coleoptera_superfamily_collapsed_tree.pdf", 
       plot = collapsed_tree_plot_image, 
       width = 16, height = 29)

ggsave("../results/FYP_figs/coleoptera_superfamily_collapsed_tree.png", 
       plot = collapsed_tree_plot_image, 
       width = 25, height = 25)
