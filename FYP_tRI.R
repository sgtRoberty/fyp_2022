library(readr)
library(dplyr)
library(ggplot2)

##### 5682-taxon aa mitogenomes #####
allmito_tri_family <- read_csv("../data/allmito_TRI_family.csv")
allmito_tri_family <- filter(allmito_tri_family,
                             ! taxon %in% c("Chrysopidae",
                                            "Coniopterygidae",
                                            "Myrmeleontidae",
                                            "Osmylidae",
                                            "Stylopidae",
                                            "Xenidae"))
allmito_tri_superfamily <- read_csv("../data/allmito_TRI_superfamily.csv")
allmito_tri_suborder <- read_csv("../data/allmito_TRI_suborder.csv")
allmito_tri_suborder <- filter(allmito_tri_suborder, 
                               ! taxon %in% c("Hemerobiiformia", 
                                              "Myrmeleontiformia",
                                              "Stylopidia"))

mean(allmito_tri_family$TRI)
mean(allmito_tri_superfamily$TRI)
mean(allmito_tri_suborder$TRI)

TRI_df <- data.frame("taxon_level" = c("Family", "Superfamily", "Suborder"),
                     "mean_tRI" = 
                         c(mean(allmito_tri_family$TRI),
                           mean(allmito_tri_superfamily$TRI),
                           mean(allmito_tri_suborder$TRI)))

meanTRI_plot <- ggplot(TRI_df, aes(x=taxon_level, y=mean_tRI, 
                                   fill=taxon_level)) +
    geom_bar(stat="identity", width=0.75, fill="#919191") + ylim(0, 1) +
    labs(x = "Taxon level", y = "Mean taxonomic retention index (tRI)") +
    geom_text(aes(label=round(mean_tRI, digits = 3)), 
              vjust=-0.6, color="black", size=5)+
    scale_x_discrete(limits=c("Family", "Superfamily", "Suborder")) + 
    theme_minimal() + 
    theme(legend.position="none",
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(vjust = 2.5),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 2.5))

ggsave("../results/FYP_figs/tRI_allmito_5682taxa_tree.pdf", 
       plot = meanTRI_plot, width = 7.5, height = 6.5)


##### RapidNJ 5682-taxon aa mitogenomes #####
RapidNJ_allmito_tri_family <- read_csv("../data/RapidNJ_noneg_TRI_family.csv")
RapidNJ_allmito_tri_family <- filter(RapidNJ_allmito_tri_family,
                                     ! taxon %in% c("Chrysopidae",
                                                    "Coniopterygidae",
                                                    "Myrmeleontidae",
                                                    "Osmylidae",
                                                    "Stylopidae",
                                                    "Xenidae"))
RapidNJ_allmito_tri_superfamily <- read_csv("../data/RapidNJ_noneg_TRI_superfamily.csv")
RapidNJ_allmito_tri_suborder <- read_csv("../data/RapidNJ_noneg_TRI_suborder.csv")
RapidNJ_allmito_tri_suborder <- filter(RapidNJ_allmito_tri_suborder, 
                                       ! taxon %in% c("Hemerobiiformia", 
                                                      "Myrmeleontiformia",
                                                      "Stylopidia"))

mean(RapidNJ_allmito_tri_family$TRI)
mean(RapidNJ_allmito_tri_superfamily$TRI)
mean(RapidNJ_allmito_tri_suborder$TRI)

RapidNJ_TRI_df <- data.frame("taxon_level" = c("Family", "Superfamily", "Suborder"),
                             "mean_tRI" = 
                                 c(mean(RapidNJ_allmito_tri_family$TRI),
                                   mean(RapidNJ_allmito_tri_superfamily$TRI),
                                   mean(RapidNJ_allmito_tri_suborder$TRI)))

RapidNJ_meanTRI_plot <- ggplot(RapidNJ_TRI_df, aes(x=taxon_level, y=mean_tRI, 
                                                   fill=taxon_level)) +
    geom_bar(stat="identity", width=0.75, fill="#919191") + ylim(0, 1) +
    labs(x = "Taxon level", y = "Mean taxonomic retention index (tRI)") +
    geom_text(aes(label=round(mean_tRI, digits = 3)), 
              vjust=-0.6, color="black", size=5)+
    scale_x_discrete(limits=c("Family", "Superfamily", "Suborder")) + 
    theme_minimal() +
    theme(legend.position="none",
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(vjust = 2.5),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 2.5))

ggsave("../results/FYP_figs/RapidNJ_tRI_allmito_5682taxa_tree.pdf", 
       plot = RapidNJ_meanTRI_plot, width = 7.5, height = 6.5)


##### Global Malaise trap data #####
allmito_malaise_tri_family <- read_csv("../data/5709taxa_malaise_tri_family.csv")
png("../results/tRI_5709taxa_globalmalaise_tree.png",
    width = 700, height = 700)
text(barplot(mean(allmito_malaise_tri_family$TRI),
             main = "Best tree of 5709 nt mitogenomes + global Malaise trap data",
             ylab = "taxonomic RI",
             xlab = "Taxon level",
             names.arg = "Family",
             ylim = c(0, 1.1)),
     y = mean(allmito_malaise_tri_family$TRI),
     labels = mean(allmito_malaise_tri_family$TRI),
     pos = 3)
dev.off()
